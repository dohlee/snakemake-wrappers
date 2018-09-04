__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

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
output = snakemake.output[0]

# Extract parameters.
# Extract optional parameters.
extra = snakemake.params.get('extra', '')
file_format_option = ''
if alignment.endswith('.bam') and '-f bam' not in extra:
    file_format_option = '-f bam'

mode_option = optionify_params('mode', '--mode')
stranded_option = optionify_params('stranded', '--stranded')

# Execute shell command.
shell(
    "("
    "htseq-count "
    "{extra} "
    "{file_format_option} "
    "{mode_option} "
    "{stranded_option} "
    "{alignment} "
    "{annotation}"
    ") > {output} "
    "{log}"
)
