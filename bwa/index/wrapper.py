__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from os import path
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

# Define exception classes.
class RuleInputException(Exception):
    pass

class RuleParameterException(Exception):
    pass

class RuleOutputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('a', '-a'))
user_parameters.append(optionify_params('b', '-b'))
user_parameters.append(optionify_params('6', '-6'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Assert input and output have been correctly given.
if len(snakemake.input) != 1:
    raise RuleInputException('Please check your reference genome has been correctly given. It should be given as a single file.')

if len(snakemake.output) != 5:
    raise RuleOutputException('bwa-mem generates 5 outputs, *.amb, *.ann, *.bwt, *.pac, and *.sa. Please check your specified output.')

# Extract required inputs.
reference = snakemake.input[0]
prefix = path.splitext(reference)[0]

# Execute shell command.
shell(
    "("
    "bwa index "
    "-p {prefix} "
    "{extra} "
    "{user_parameters} "
    "{reference} "
    ") "
    "{log}"
)
