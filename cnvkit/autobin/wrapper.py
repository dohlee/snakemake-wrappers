__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import os

from snakemake.shell import shell

def is_defined_by_user(*params):
    extra = snakemake.params.get('extra', '')
    for param in params:
        if param in extra:
            return True
    return False

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + str(snakemake.params[parameter])
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
bams = snakemake.input.bams
targets = snakemake.input.targets
access = snakemake.input.access

# Extract optional parameters.
user_parameters = []
user_parameters = ' '.join(user_parameters)

# Execute shell command.
shell(
    "("
    "cnvkit.py autobin {bams} "
    "--targets {targets} "
    "--access {access} "
    "{extra} "
    ") "
    "{log}"
)
