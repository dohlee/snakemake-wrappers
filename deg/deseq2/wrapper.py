__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import rpy2.robjects as ro
from rpy2.robjects import r
import pandas as pd
import logging
import warnings

from os import path
from snakemake.shell import shell
from rpy2.robjects import pandas2ri
from rpy2.robjects import Formula
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
from rpy2.rinterface import RRuntimeWarning

warnings.filterwarnings('ignore', category=RRuntimeWarning)
base = importr('base')

# Logger settings.
LOGGER_NAME = 'DESeq2'
logger = logging.getLogger(LOGGER_NAME)
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()
formatter = logging.Formatter('[%(levelname)-.1s %(asctime)s %(name)s] %(message)s')
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


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
    logging.debug('Testing convergence...%g, %g' % (vec[l], vec[l-1]))
    return abs(vec[l] - vec[l-1]) < threshold


# Load required package.
deseq2 = import_bioc('DESeq2')
biocparallel = import_bioc('biocparallel')
biocparallel.register(MulticoreParam(snakemake.threads))

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
data = pd.read_table(snakemake.input.data, index_col=0)  # Input Gene-by-Sample raw count data.
condition = pd.read_table(snakemake.input.condition, index_col=0, names=['condition'])  # Input condition file which indicates to which condition each sample belongs.
logger.info('%d(genes) x %d(samples) data matrix and %d sample conditions are given.' % (data.shape[0], data.shape[1], len(condition.index)))
logger.info('Headers: %s...' % ' '.join(data.columns[:3]))
logger.info('Gene identifiers: %s...' % ' '.join(data.index[:3]))

intersecting_samples = [sample for sample in data.columns if sample in condition.index]
data = data[intersecting_samples]
condition = condition.loc[intersecting_samples]
logger.info('%d samples will be used for DEG discovery.' % len(intersecting_samples))

pandas2ri.activate()
r_data_matrix = r['data.matrix'](pandas2ri.py2ri(data))
pandas2ri.deactivate()

logger.info('Making DESeq dataset from matrix.')
r_deseq_dataset = deseq2.DESeqDataSetFromMatrix(countData=r_data_matrix, colData=condition, design=Formula('~ condition'))

logger.info('Discovering DEGs.')
logger.info('Running DESeq.')
r_deseq_result = deseq2.DESeq(r_deseq_dataset)
r_overall_result = deseq2.results(r_deseq_result)

# Write DEGs to output file.
pandas2ri.activate()
result = pd.DataFrame(pandas2ri.ri2py(r_overall_result))
pandas2ri.deactivate()

result.to_csv(snakemake.output.result, sep='\t')
