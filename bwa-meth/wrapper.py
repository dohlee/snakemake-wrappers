__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

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
user_parameters.append(optionify_params('read_group', '--read-group'))
user_parameters.append(optionify_params('set_as_failed', '--set-as-failed'))
user_parameters.append(optionify_params('interleaved', '--interleaved'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
reads = snakemake.input.reads
reference = snakemake.input.reference

# Extract required outputs.
output = snakemake.output[0]

pipe_command = ''
# bwa-mem output defaults to sam. Convert sam to bam with samtools.
# TODO: Use sambamba with appropriate number of threads.
if output.endswith('.bam'):
    pipe_command = '| samtools view -Sb -'

# Execute shell command.
shell(
    "("
    "bwameth.py "
    "{extra} "
    "{user_parameters} "
    "-t {snakemake.threads} "
    "--reference {reference} "
    "{reads} "
    "{pipe_command} > "
    "{output} "
    ")"
    "{log}"
)
