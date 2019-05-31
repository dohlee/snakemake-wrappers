__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
from snakemake.shell import shell

# Define utility function.
def get_common_prefixes(strings):
    all_same = lambda x: all(x[0] == y for y in x)

    prefix_tuples = itertools.takewhile(all_same, zip(*strings))
    return ''.join(x[0] for x in prefix_tuples).strip('.')

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        if str(snakemake.params[parameter]) == '':
            return ''
        if type(snakemake.params[parameter]) == bool:
            return option
        else:
            return option + ' ' + str(snakemake.params[parameter])
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
input_file = snakemake.input[0]
input_command = '-i %s' % input_file

# Extract optional inputs.
output_file = snakemake.output[0]
output_command = '-o %s' % output_file

# Extract user parameters.
user_parameters = []
user_parameters.append(optionify_params('gsize', '--gsize'))
user_parameters.append(optionify_params('tsize', '--tsize'))
user_parameters.append(optionify_params('pvalue', '--pvalue'))
user_parameters.append(optionify_params('keep_dup', '--keep_dup'))
user_parameters.append(optionify_params('verbose', '--verbose'))
user_parameters = ' '.join([p for p in user_parameters if not p != ''])
    
# Execute shell command.
shell(
    "("
    "macs2 filterdup "
    "{input_command} "
    "{output_command} "
    "{user_parameters} "
    "{extra} "
    ") "
    "{log}"
)
