__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


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
LOGGER_NAME = 'edgeR'
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
dispersion = snakemake.params.get('dispersion', 'common')  # common, trended, tagwise
assert dispersion in ['common', 'trended', 'tagwise'], 'Dispersion estimation should be one of ["common", "trended", "tagwise"].'
if verbose:
    logger.setLevel(cleanlog.DEBUG)
else:
    logger.setLevel(cleanlog.INFO)

# Load required package.
edgeR = import_bioc('edgeR')

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

logger.debug('Making DGE List.')
dge_list = edgeR.DGEList(counts=data, genes=r['rownames'](data))

logger.debug('Filtering genes that are not expressed in many of the samples.')
counts_per_million = edgeR.cpm(dge_list)
count_check = counts_per_million.ro > 1
keep = r['which'](r['rowSums'](count_check).ro >= 2)
dge_list = dge_list.rx(keep, True)

logger.debug('Normalizing data.')
dge_list = r['calcNormFactors'](dge_list, method="TMM")

logger.debug('Setting up the model.')
formula = Formula('~ condition$condition')
env = formula.environment
env['condition'] = condition
design_matrix = r['model.matrix'](formula)

logger.debug('Estimating dispersions.')
if dispersion == 'common':
    dge_list = edgeR.estimateGLMCommonDisp(dge_list, design=design_matrix)
elif dispersion == 'trended':
    dge_list = edgeR.estimateGLMTrendedDisp(dge_list, design=design_matrix)
else:
    dge_list = edgeR.estimateGLMTagwiseDisp(dge_list, design=design_matrix)

logger.debug('Finding DEGs.')
fit = edgeR.glmFit(dge_list, design_matrix)
lrt = edgeR.glmLRT(fit)

edger_result = edgeR.topTags(lrt, n=len(data.index))
edger_result_table = pandas2ri.ri2py(edger_result[edger_result.names.index('table')])
edger_result_table.to_csv(snakemake.output.result, sep='\t')

with open(snakemake.output.deg_list, 'w') as outFile:
    print('\n'.join(edger_result_table[edger_result_table.PValue < cutoff].genes.values), file=outFile)
