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
output = snakemake.output

# Extract optional parameters.
annotation_format = optionify_params('annotation_format', '-F')
feature_type = optionify_params('feature_type', '-t')
attribute_type = optionify_params('attribute_type', '-g')
min_overlap = optionify_params('min_overlap', '--minOverlap')

# Execute shell command.
shell(
    "("
    "featureCounts "
    "{extra} "
    "{annotation_format} "
    "{attribute_type} "
    "{min_overlap} "
    "-T {snakemake.threads} "
    "-a {annotation} "
    "-o {output} "
    "{alignment} "
    ") "
    "{log}"
)
