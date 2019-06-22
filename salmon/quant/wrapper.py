__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
from snakemake.shell import shell

# Define utility function.
def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        if str(snakemake.params[parameter]) == '':
            return ''
        if type(snakemake.params[parameter]) == bool:
            if snakemake.params[parameter]:
                return option
            else:
                return ''
        else:
            return option + ' ' + str(snakemake.params[parameter])
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
library_type_parameter = optionify_params('libType', '--libType')
user_parameters = []
user_parameters.append(optionify_params('seqBias', '--seqBias'))
user_parameters.append(optionify_params('gcBias', '--gcBias'))
user_parameters.append(optionify_params('incompatPrior', '--incompatPrior'))
user_parameters.append(optionify_params('geneMap', '--geneMap'))
user_parameters.append(optionify_params('meta', '--meta'))
user_parameters.append(optionify_params('discardOrphansQuasi', '--discardOrphansQuasi'))
user_parameters.append(optionify_params('validateMappings', '--validateMappings'))
user_parameters.append(optionify_params('consensusSlack', '--consensusSlack'))
user_parameters.append(optionify_params('minScoreFraction', '--minScoreFraction'))
user_parameters.append(optionify_params('maxMMPExtension', '--maxMMPExtension'))
user_parameters.append(optionify_params('ma', '--ma'))
user_parameters.append(optionify_params('mp', '--mp'))
user_parameters.append(optionify_params('go', '--go'))
user_parameters.append(optionify_params('ge', '--ge'))
user_parameters.append(optionify_params('bandwidth', '--bandwidth'))
user_parameters.append(optionify_params('allowDovetail', '--allowDovetail'))
user_parameters.append(optionify_params('recoverOrphans', '--recoverOrphans'))
user_parameters.append(optionify_params('mimicBT2', '--mimicBT2'))
user_parameters.append(optionify_params('mimicStrictBT2', '--mimicStrictBT2'))
user_parameters.append(optionify_params('hardFilter', '--hardFilter'))
user_parameters.append(optionify_params('writeMappings', '--writeMappings'))
user_parameters.append(optionify_params('consistentHits', '--consistentHits'))
user_parameters.append(optionify_params('numBootstraps', '--numBootstraps'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
reads = snakemake.input.reads
index = snakemake.input.index

if len(reads) == 1: # Single-end
    read_command = '-r %s' % reads[0]
elif len(reads) == 2:  # Paired-end
    read_command = '-1 %s -2 %s' % (reads[0], reads[1])
else:
    raise ValueError('Please give one or two reads.')

# Extract required outputs.
quant = snakemake.output.quant
lib = snakemake.output.lib
output_directory = path.dirname(quant)

# Execute shell command.
shell(
    "("
    "salmon quant "
    "{library_type_parameter} "
    "{read_command} "
    "-i {index} "
    "-o {output_directory} "
    "-p {snakemake.threads} "
    "{extra} "
    "{user_parameters} "
    ")"
    "{log}"
)
