__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import rpy2.robjects as ro
from rpy2.robjects import r
import pandas as pd
import cleanlog
import warnings

from os import path
from snakemake.shell import shell
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
from rpy2.rinterface import RRuntimeWarning

warnings.filterwarnings('ignore', category=RRuntimeWarning)
base = importr('base')
pandas2ri.activate()

# Logger settings.
LOGGER_NAME = 'EBSeq'
logger = cleanlog.ColoredLogger(LOGGER_NAME)


def import_bioc(package_name):
    """Import bioconductor packages. If not installed, install the package."""
    base.source("http://www.bioconductor.org/biocLite.R")
    biocinstaller = importr("BiocInstaller")
    try:
        logger.info('Importing %s.' % package_name)
        package = importr(package_name)
    except RRuntimeError:
        logger.info('%s is not installed. Installing.' % package_name)
        biocinstaller.biocLite(package_name)
        logger.info('Importing %s again.' % package_name)
        package = importr(package_name)

    logger.info('Package %s has been successfully imported.' % package_name)
    return package


def converged(vec, l, threshold=1e-3):
    logger.debug('Testing convergence...%g, %g' % (vec[l], vec[l-1]))
    return abs(vec[l] - vec[l-1]) < threshold

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
cutoff = snakemake.params.get('cutoff', 0.05)
verbose = snakemake.params.get('verbose', False)
if verbose:
    logger.setLevel(cleanlog.DEBUG)
else:
    logger.setLevel(cleanlog.INFO)

# Load required package.
ebseq = import_bioc('EBSeq')

# Extract required arguments.
data = pd.read_table(snakemake.input.data, index_col=0)  # Input Gene-by-Sample raw count data.
condition = pd.read_table(snakemake.input.condition, index_col=0, names=['condition'])  # Input condition file which indicates to which condition each sample belongs.
logger.info('%d(genes) x %d(samples) data matrix and %d sample conditions are given.' % (data.shape[0], data.shape[1], len(condition.index)))
logger.debug('Headers: %s...' % ' '.join(data.columns[:3]))
logger.debug('Gene identifiers: %s...' % ' '.join(data.index[:3]))

intersecting_samples = [sample for sample in data.columns if sample in condition.index]
data = data[intersecting_samples]

condition = list(condition.loc[intersecting_samples].condition.values)
logger.info('%d samples will be used for DEG discovery.' % len(intersecting_samples))


r_data_matrix = r['data.matrix'](pandas2ri.py2ri(data))
r_samples = r.colnames(r_data_matrix)
r_conditions = ro.FactorVector(condition)

logger.debug('Computing size factors.')
r_size_factors = ebseq.MedianNorm(r_data_matrix)

logger.info('Discovering DEGs.')
logger.info('Running EBTest.')

num_iteration = 0
while True:
    # Increase iteration numbers if the conditons are not met.
    # Hopefully most of the tie, 10 iterations will be enough for convergence.
    num_iteration += 10
    r_eb_out = ebseq.EBTest(Data=r_data_matrix, Conditions=r_conditions, sizeFactors=r_size_factors, maxround=num_iteration)
    logger.info('Running GetDEResults. (FDR cutoff = %.3f)' % cutoff)
    r_eb_de_result = ebseq.GetDEResults(r_eb_out, FDR=cutoff)

    # Check convergences.
    # Each parameter should change less than 1e-3 between the last two iterations.
    l = len(r_eb_out[r_eb_out.names.index('Alpha')]) - 1
    if converged(r_eb_out[r_eb_out.names.index('Alpha')], l) and converged(r_eb_out[r_eb_out.names.index('Beta')], l) and converged(r_eb_out[r_eb_out.names.index('P')], l):
        break

# Ouput DEG list.
r_deg_found = r_eb_de_result[r_eb_de_result.names.index('DEfound')]
deg_found = pd.DataFrame(pandas2ri.ri2py(r_deg_found))
logger.info('Writing ebseq DEG list to %s' % snakemake.output.deg_list)
deg_found.to_csv(snakemake.output.deg_list, sep='\t', index=False, header=None)

# Output fold-chanages.
r_fold_changes = ebseq.PostFC(r_eb_out)
genes = list(r_fold_changes[r_fold_changes.names.index('PostFC')].names)
logger.debug('genes: %d' % len(genes))
post_fc = list(r_fold_changes[r_fold_changes.names.index('PostFC')])
logger.debug('post_fc: %d' % len(post_fc))
real_fc = list(r_fold_changes[r_fold_changes.names.index('RealFC')])
logger.debug('real_fc: %d' % len(real_fc))
direction = list(r_fold_changes[r_fold_changes.names.index('Direction')])
fold_changes = pd.DataFrame({'PostFC':post_fc, 'RealFC':real_fc, 'Direction':direction}, index=genes)
fold_changes.index.name = 'Gene ID'

logger.info('Writing ebseq fold change to %s' % snakemake.output.fold_change)
fold_changes.to_csv(snakemake.output.fold_change, sep='\t')

# Output result table.
r_result = r_eb_de_result[r_eb_de_result.names.index('PPMat')]
logger.debug('r_result: ', type(pandas2ri.ri2py(r_result)), pandas2ri.ri2py(r_result).shape)
result = pd.DataFrame(pandas2ri.ri2py(r_result), index=list(r['rownames'](r_result)), columns=['padj', '1-padj'])
result.index.name = 'Gene ID'

logger.info('Writing ebseq result to %s.' % snakemake.output.result)
result.to_csv(snakemake.output.result, sep='\t')
