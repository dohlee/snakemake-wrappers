__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

class RuleOutputException(Exception):
    pass

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
bam = snakemake.input.bam
reference = snakemake.input.reference

# Extract required outputs.
output = snakemake.output[0]

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params('mapping_quality_threshold', '-q'))
user_parameters.append(optionify_params('sequencing_quality_threshold', '-p'))
user_parameters = ' '.join(user_parameters)

wrapper_parameters = []
if not is_defined_by_user('--noSVG'):
    wrapper_parameters.append('--noSVG')
wrapper_parameters = ' '.join(wrapper_parameters)

pipe_command = '> %s' % output

# Execute shell command.
shell(
    "("
    "MethylDackel mbias "
    "-@ {snakemake.threads} "
    "{extra} "
    "{user_parameters} "
    "{wrapper_parameters} "
    "{reference} "
    "{bam} "
    "{pipe_command} "
    ") "
    "{log}"
)
