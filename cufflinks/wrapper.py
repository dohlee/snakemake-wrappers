__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.input[parameter])
        if param:
            return option + ' ' + str(snakemake.input[parameter])
        else:
            return ''
    except AttributeError:
        return ''

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.params[parameter])
        if param:
            return option + ' ' + str(snakemake.params[parameter])
        else:
            return ''
    except AttributeError:
        return ''

# Extract required inputs.
alignment = snakemake.input.alignment
annotation = snakemake.input.annotation

# Extract required outputs.
output = snakemake.output
output_directory = os.path.dirname(output[0])

# Extract optional parameters.
annotation_format = optionify_params('annotation_format', '-F')
random_seed = optionify_params('random_seed', '--seed')

# Execute shell command.
shell(
    "("
    "cufflinks "
    "{extra} "
    "{random_seed} "
    "-p {snakemake.threads} "
    "-G {annotation} "
    "-o {output_directory} "
    "{alignment} "
    ") "
    "{log}"
)
