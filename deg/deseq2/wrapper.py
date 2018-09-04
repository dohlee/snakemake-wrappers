__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import rpy2.robjects as ro
from rpy2.robjects import r
import pandas as pd
import cleanlog
import warnings

from rpy2.robjects import pandas2ri
from rpy2.robjects import Formula
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
from rpy2.rinterface import RRuntimeWarning

warnings.filterwarnings('ignore', category=RRuntimeWarning)
base = importr('base')
pandas2ri.activate()

# Logger settings.
LOGGER_NAME = 'DESeq2'
logger = cleanlog.ColoredLogger(LOGGER_NAME)


def import_bioc(package_name):
    """Import bioconductor packages. If not installed, install the package."""
    base.source("http://www.bioconductor.org/biocLite.R")
    biocinstaller = importr("BiocInstaller")
    try:
        logger.debug('Importing %s.' % package_name)
        package = importr(package_name)
    except RRuntimeError:
        logger.warning('%s is not installed. Installing.' % package_name)
        biocinstaller.biocLite(package_name)
        logger.debug('Importing %s again.' % package_name)
        package = importr(package_name)

    logger.debug('Package %s has been successfully imported.' % package_name)
    return package

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
deseq2 = import_bioc('DESeq2')
biocparallel = import_bioc('BiocParallel')
# Use threading to speed things up.
biocparallel.register(biocparallel.MulticoreParam(snakemake.threads))

# Extract required arguments.
data = pd.read_table(snakemake.input.data, index_col=0)  # Input Gene-by-Sample raw count data.
condition = pd.read_table(snakemake.input.condition, index_col=0, names=['condition'])  # Input condition file which indicates to which condition each sample belongs.
logger.info('%d(genes) x %d(samples) data matrix and %d sample conditions are given.' % (data.shape[0], data.shape[1], len(condition.index)))
logger.debug('Headers: %s...' % ' '.join(data.columns[:3]))
logger.debug('Gene identifiers: %s...' % ' '.join(data.index[:3]))

intersecting_samples = [sample for sample in data.columns if sample in condition.index]
data = data[intersecting_samples]
condition = condition.loc[intersecting_samples]
logger.info('%d samples will be used for DEG discovery.' % len(intersecting_samples))

# NOTE: Converting pandas dataframe into R data.matrix seems to convert hyphens into dots.
# Since sample names in countData and colData should match, convert hyphens in sample names of condition into dots.
condition.index = [sample.replace('-', '.') for sample in condition.index]

logger.debug('Making DESeq dataset from matrix.')
r_deseq_dataset = deseq2.DESeqDataSetFromMatrix(countData=data, colData=condition, design=Formula('~ condition'))

logger.debug('Discovering DEGs.')
logger.info('Running DESeq2 with %d threads.' % snakemake.threads)
r_deseq_result = deseq2.DESeq(r_deseq_dataset)
r_overall_result = r['as.data.frame'](deseq2.results(r_deseq_result))
# Write DEGs to output file.
result = pd.DataFrame(pandas2ri.ri2py(r_overall_result))
result.index = list(pandas2ri.ri2py(r['rownames'](r_overall_result)))
result.index.name = 'Gene ID'

logger.info('Writing deseq2 result to %s.' % snakemake.output.result)
result.to_csv(snakemake.output.result, sep='\t')

# Export found DEGs into a list.
logger.info('Writing deseq2 DEG list(padj cutoff=%.3f) to %s.' % (cutoff, snakemake.output.deg_list))
degs = result[result.padj < cutoff]
with open(snakemake.output.deg_list, 'w') as outFile:
    print('\n'.join([s for s in degs.index]), file=outFile)
