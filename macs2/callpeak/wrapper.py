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

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
treatment = snakemake.input.treatment
input_command = '-t %s' % treatment

# Extract optional inputs.
control = snakemake.input.control
if control:
    input_command += ' -c %s' % control

# Extract required outputs.
output_dir = path.dirname(snakemake.output.peak)
output_command = '--outdir %s' % output_dir

# Extract common prefixes of output files.
common_prefix = get_common_prefixes([path.basename(f) for f in snakemake.output])[:-1]
name_command = '-n %s' % common_prefix

# Extract user parameters.
user_parameters = []
user_parameters.append(optionify_params('genome_size', '-g'))
user_parameters.append(optionify_params('broad', '--broad'))
user_parameters.append(optionify_params('seed', '--seed'))
user_parameters.append(optionify_params('bedGraph_out', '--bdg'))
user_parameters.append(optionify_params('q_value_cutoff', '-q'))
user_parameters.append(optionify_params('p_value_cutoff', '-p'))
user_parameters = ' '.join(user_parameters)
    
# Execute shell command.
shell(
    "("
    "macs2 callpeak "
    "{input_command} "
    "{output_command} "
    "{name_command} "
    "{user_parameters} "
    "{extra} "
    ") "
    "{log}"
)
