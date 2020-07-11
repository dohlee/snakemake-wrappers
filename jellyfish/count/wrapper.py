__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2020, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path
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
user_parameters.append(optionify_params('mer_len', '--mer-len'))
user_parameters.append(optionify_params('size', '--size'))
user_parameters.append(optionify_params('sam', '--sam'))
user_parameters.append(optionify_params('Files', '--Files'))
user_parameters.append(optionify_params('generator', '--generator'))
user_parameters.append(optionify_params('Generators', '--Generators'))
user_parameters.append(optionify_params('shell_', '--shell'))
user_parameters.append(optionify_params('counter_len', '--counter-len'))
user_parameters.append(optionify_params('out_counter_len', '--out-counter-len'))
user_parameters.append(optionify_params('canonical', '--canonical'))
user_parameters.append(optionify_params('bc', '--bc'))
user_parameters.append(optionify_params('bf_size', '--bf-size'))
user_parameters.append(optionify_params('bf_fp', '--bf-fp'))
user_parameters.append(optionify_params('if_', '--if'))
user_parameters.append(optionify_params('min_qual_char', '--min-qual-char-string'))
user_parameters.append(optionify_params('quality_start', '--quality-start'))
user_parameters.append(optionify_params('min_quality', '--min-quality'))
user_parameters.append(optionify_params('reprobes', '--reprobes'))
user_parameters.append(optionify_params('text', '--text'))
user_parameters.append(optionify_params('dist', '--dist'))
user_parameters.append(optionify_params('lower_count', '--lower-count'))
user_parameters.append(optionify_params('upper_count', '--upper-count'))
user_parameters.append(optionify_params('timing', '--timing'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
fastq = snakemake.input.fastq

# Compose input command.
if fastq[0].endswith('.gz'):
    input_command = 'gunzip -c ' + ' '.join(list(snakemake.input.fastq))
else:
    input_command = 'cat ' + ' '.join(list(snakemake.input.fastq))

# Execute shell command.
shell(
    "("
    "jellyfish count "
    "{extra} "
    "{user_parameters} "
    "-t {snakemake.threads} "
    "-o {snakemake.output} "
    "<({input_command}) "
    ")"
    "{log}"
)
