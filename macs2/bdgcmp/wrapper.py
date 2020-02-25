__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2020, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
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

def determine_methods(outputs):
    lookup = {
        '_ppois.bdg': 'ppois',
        '_qpois.bdg': 'qpois',
        '_subtract.bdg': 'subtract',
        '_FE.bdg': 'FE',
        '_logFE.bdg': 'logFE',
        '_logLR.bdg': 'logLR',
        '_slogLR.bdg': 'slogLR',
    }

    methods = []
    for output in outputs:
        found_suffix = False

        for suffix, method in lookup.items():
            if output.endswith(suffix):
                methods.append(method) 
                found_suffix = True

        if not found_suffix:
            raise ValueError('Invalid output file name: ' + output)

    return '--method ' + ' '.join(methods)

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
treatment = snakemake.input.treatment
control = snakemake.input.control
input_command = '-t %s -c %s' % (treatment, control)

# Extract required outputs.
output = list(snakemake.output)
output_command = '-o %s' % ' '.join(output)

# Extract user parameters.
user_parameters = []
user_parameters.append(optionify_params('scaling_factor', '--scaling-factor'))
user_parameters.append(optionify_params('pseudocount', '--pseudocount'))
user_parameters = ' '.join(user_parameters)

method_parameters = determine_methods(output)
    
# Execute shell command.
shell(
    "("
    "macs2 bdgcmp "
    "{input_command} "
    "{output_command} "
    "{user_parameters} "
    "{method_parameters} "
    "{extra} "
    ") "
    "{log}"
)
