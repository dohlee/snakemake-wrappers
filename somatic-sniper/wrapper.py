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

def optionify_params(parameter, option, default=None):
    """Return optionified parameter."""
    try:
        return option + ' ' + snakemake.params[parameter]
    except AttributeError:
        return '' if default is None else option + ' ' + str(default)

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
tumor_bam = snakemake.input.tumor_bam
normal_bam = snakemake.input.normal_bam
reference = snakemake.input.reference

# Extract required outputs.
output = snakemake.output[0]

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params(mapping_quality_cutoff, '-q', default=20))
user_parameters.append(optionify_params(calling_quality_cutoff, '-Q', default=15))
user_parameters.append(optionify_params(tumor_sample_name, '-t'))
user_parameters.append(optionify_params(normal_sample_name, '-n'))
user_parameters = ' '.join(user_parameters)

wrapper_parameters = []
if not is_defined_by_user('-F'):
    wrapper_parameters.append('-F vcf')
wrapper_parameters = ' '.join(wrapper_parameters)

# Execute shell command.
shell(
    "("
    "bam-somaticsniper "
    "{extra} "
    "{user_parameters} "
    "{wrapper_parameters} "
    "-f {reference} "
    "{tumor_bam} "
    "{normal_bam} "
    "{output} "
    ") "
    "{log}"
)
