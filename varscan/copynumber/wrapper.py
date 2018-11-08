__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


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

# Extract logself.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
tumor_bam = snakemake.input.tumor_bam
normal_bam = snakemake.input.normal_bam
reference = snakemake.input.reference
output_prefix = snakemake.params.output_prefix

q_cutoff = snakemake.params.get('pileup_quality_cutoff', 1)
pipe_command = f'samtools mpileup -q {q_cutoff} -f {reference} {normal_bam} {tumor_bam} |'

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params('min_base_qual', '--min-base-qual'))
user_parameters.append(optionify_params('min_map_qual', '--min-map-qual'))
user_parameters.append(optionify_params('min_coverage', '--min-coverage'))
user_parameters.append(optionify_params('min_segment_size', '--min-segment-size'))
user_parameters.append(optionify_params('max_segment_size', '--max-segment-size'))
user_parameters.append(optionify_params('p_value', '--p-value'))
user_parameters.append(optionify_params('data_ratio', '--data-ratio'))
user_parameters = ' '.join(user_parameters)

# Execute shell command.
shell(
    "("
    "{pipe_command} "
    "varscan copynumber "
    "{output_prefix} "
    "--mpileup 1 "
    "{extra} "
    "{user_parameters} "
    ") "
    "{log}"
)
