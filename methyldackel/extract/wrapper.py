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
user_parameters = []
user_parameters.append(optionify_params('q', '-q'))
user_parameters.append(optionify_params('p', '-p'))
user_parameters.append(optionify_params('D', '-D'))
user_parameters.append(optionify_params('d', '-d'))
user_parameters.append(optionify_params('r', '-r'))
user_parameters.append(optionify_params('l', '-l'))
user_parameters.append(optionify_params('keepStrand', '--keepStrand'))
user_parameters.append(optionify_params('chunkSize', '--chunkSize'))
user_parameters.append(optionify_params('mergeContext', '--mergeContext'))
user_parameters.append(optionify_params('opref', '--opref'))
user_parameters.append(optionify_params('keepDupes', '--keepDupes'))
user_parameters.append(optionify_params('keepSingleton', '--keepSingleton'))
user_parameters.append(optionify_params('keepDiscordant', '--keepDiscordant'))
user_parameters.append(optionify_params('ignoreFlags', '--ignoreFlags'))
user_parameters.append(optionify_params('requireFlags', '--requireFlags'))
user_parameters.append(optionify_params('noCpG', '--noCpG'))
user_parameters.append(optionify_params('CHG', '--CHG'))
user_parameters.append(optionify_params('CHH', '--CHH'))
user_parameters.append(optionify_params('minOppositeDepth', '--minOppositeDepth'))
user_parameters.append(optionify_params('maxVariantFrac', '--maxVariantFrac'))
user_parameters.append(optionify_params('methylKit', '--methylKit'))
user_parameters.append(optionify_params('cytosine_report', '--cytosine_report'))
user_parameters.append(optionify_params('OT', '--OT'))
user_parameters.append(optionify_params('OB', '--OB'))
user_parameters.append(optionify_params('CTOT', '--CTOT'))
user_parameters.append(optionify_params('CTOB', '--CTOB'))
user_parameters.append(optionify_params('nOT', '--nOT'))
user_parameters.append(optionify_params('nOB', '--nOB'))
user_parameters.append(optionify_params('nCTOT', '--nCTOT'))
user_parameters.append(optionify_params('nCTOB', '--nCTOB'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
bam = snakemake.input.bam
reference = snakemake.input.reference

# Extract required outputs.
output = snakemake.output[0]
if output.endswith('_CpG.bedGraph'):
    type_option = ''
elif output.endswith('_CpG.meth.bedGraph'):
    type_option = '--fraction'
elif output.endswith('_CpG.counts.bedGraph'):
    type_option = '--counts'
elif output.endswith('_CpG.logit.bedGraph'):
    type_option = '--logit'
else:
    raise ValueError('Unrecognized output format: %s. Use *_CpG.bedGrah or *_CpG.{meth,counts,logit}.bedGraph')

# Execute shell command.
shell(
    "("
    "MethylDackel extract "
    "-@ {snakemake.threads} "
    "{extra} "
    "{user_parameters} "
    "{type_option} "
    "{reference} "
    "{bam} "
    ")"
    "{log}"
)
