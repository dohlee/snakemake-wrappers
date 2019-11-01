__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import sys
from os import path

from snakemake.shell import shell

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + snakemake.input[parameter]
    except AttributeError:
        return ''

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

# Define exception classes.
class RuleParameterException(Exception):
    pass

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
java_options = snakemake.params.get('java_options', '-Xmx4g')

# Extract required arguments.
vcf = snakemake.input.vcf
db = snakemake.input.db

output = snakemake.output[0]
pipe_command = '> %s' % output

user_parameters = []
user_parameters.append(optionify_params('id', '-id'))
user_parameters.append(optionify_params('info', '-info'))
user_parameters.append(optionify_params('tabix', '-tabix'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Execute shell command.
shell(
    "("
    "SnpSift annotate "
    "{java_options} "
    "{extra} "
    "{user_parameters} "
    "{db} "
    "{vcf} "
    "{pipe_command}"
    ")"
    "{log}"
)
