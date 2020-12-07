__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2020, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


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

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params('annotation_format', '-F'))
user_parameters.append(optionify_params('feature_type', '-t'))
user_parameters.append(optionify_params('attribute_type', '-g'))
user_parameters.append(optionify_params('extra_attributes', '--extraAttributes'))
user_parameters.append(optionify_params('chromosome_name_alias', '-A'))
user_parameters.append(optionify_params('feature_level_counting', '-f'))
user_parameters.append(optionify_params('multi_overlapping_reads', '-O'))
user_parameters.append(optionify_params('min_overlap', '--minOverlap'))
user_parameters.append(optionify_params('frac_overlap_feature', '--fracOverlapFeature'))
user_parameters.append(optionify_params('largest_overlap', '--largestOverlap'))
user_parameters.append(optionify_params('non_overlap', '--nonOverlap'))
user_parameters.append(optionify_params('non_overlap_feature', '--nonOverlapFeature'))
user_parameters.append(optionify_params('read_extension_5', '--readExtension5'))
user_parameters.append(optionify_params('read_extension_3', '--readExtension3'))
user_parameters.append(optionify_params('read2pos', '--read2pos'))
user_parameters.append(optionify_params('count_multi_mapping', '-M'))
user_parameters = ' '.join(user_parameters)

# Execute shell command.
shell(
    "("
    "featureCounts "
    "{extra} "
    "{user_parameters} "
    "-T {snakemake.threads} "
    "-a {annotation} "
    "-o {output} "
    "{alignment} "
    ") "
    "{log}"
)
