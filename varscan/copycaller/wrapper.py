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
raw_copynumber_calls = snakemake.input.raw_copynumber_calls

# Extract required outputs.
copynumber_calls = snakemake.output.copynumber_calls
homdel = snakemake.output.homdel

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params('min_coverage', '--min-coverage'))
user_parameters.append(optionify_params('min_tumor_coverage', '--min-tumor-coverage'))
user_parameters.append(optionify_params('max_homdel_coverage', '--max-homdel-coverage'))
user_parameters.append(optionify_params('amp_threshold', '--amp-threshold'))
user_parameters.append(optionify_params('del_threshold', '--del-threshold'))
user_parameters.append(optionify_params('min_region_size', '--min-region-size'))
user_parameters.append(optionify_params('recenter_up', '--recenter-up'))
user_parameters.append(optionify_params('recenter_down', '--recenter-down'))
user_parameters = ' '.join(user_parameters)

# Execute shell command.
shell(
    "("
    "varscan copyCaller "
    "{raw_copynumber_calls} "
    "--output-file {copynumber_calls} "
    "--output-homdel-file {homdel} "
    "{extra} "
    "{user_parameters} "
    ") "
    "{log}"
)
