__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from os import dirname
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

# Extract required inputs.
bed = snakemake.input.bed
genome = snakemake.input.genome
motif = snakemake.input.motif

# Extract required outputs.
homer_result = snakemake.output.homer_result
outdir = dirname(homer_result)

# Execute shell command.
shell(
    "("
    "findMotifsGenome.pl "
    "{bed} "
    "{genome} "
    "{outdir} "
    "-find {motif} "
    "-p {snakemake.threads} > {homer_result}"
    ")"
    "{log}"
)
