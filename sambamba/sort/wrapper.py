__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
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
user_parameters.append(optionify_params('memory_limit', '--memory-limit'))
user_parameters.append(optionify_params('tmpdir', '--tmpdir'))
user_parameters.append(optionify_params('sort_by_name', '--sort-by-name'))
user_parameters.append(optionify_params('compression_level', '--compression-level'))
user_parameters.append(optionify_params('uncompressed_chunks', '--uncompressed-chunks'))
user_parameters.append(optionify_params('show_progress', '--show-progress'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
bam = snakemake.input[0]

# Extract required outputs.
sorted_bam = snakemake.output[0]

# Execute shell command.
shell(
    "("
    "sambamba sort "
    "{extra} "
    "{user_parameters} "
    "-t {snakemake.threads} "
    "-o {sorted_bam} "
    "{bam} "
    ")"
    "{log}"
)
