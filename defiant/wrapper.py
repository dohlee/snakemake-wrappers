__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

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

def format_option_for_defiant_output(parameter):
    try:
        if str(snakemake.params[parameter]) == '':
            return ''
        if type(snakemake.params[parameter]) == bool:
            if snakemake.params[parameter]:
                return parameter
            else:
                return ''
        else:
            return parameter + str(snakemake.params[parameter])
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('a', '-a'))
user_parameters.append(optionify_params('b', '-b'))
user_parameters.append(optionify_params('c', '-c'))
user_parameters.append(optionify_params('CpN', '-CpN'))
user_parameters.append(optionify_params('d', '-d'))
user_parameters.append(optionify_params('debug', '-debug'))
user_parameters.append(optionify_params('D', '-D'))
user_parameters.append(optionify_params('E', '-E'))
user_parameters.append(optionify_params('f', '-f'))
user_parameters.append(optionify_params('fdr', '-fdr'))
user_parameters.append(optionify_params('G', '-G'))
user_parameters.append(optionify_params('l', '-l'))
user_parameters.append(optionify_params('N', '-N'))
user_parameters.append(optionify_params('p', '-p'))
user_parameters.append(optionify_params('P', '-P'))
user_parameters.append(optionify_params('q', '-q'))
user_parameters.append(optionify_params('r', '-r'))
user_parameters.append(optionify_params('R', '-R'))
user_parameters.append(optionify_params('s', '-s'))
user_parameters.append(optionify_params('S', '-S'))
user_parameters.append(optionify_params('U', '-U'))
user_parameters.append(optionify_params('v', '-v'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
a_files = ','.join(snakemake.input.a)
b_files = ','.join(snakemake.input.b)

# Extract required outputs.
out = snakemake.output

# Predict raw Defiant output and rename it to the specified output.
# Defiant shows the values of -c, -CpN, -d, -G, -p, -P, -S options
# in the output file name by default.
predicted_out = [f'{snakemake.wildcards.a}_vs_{snakemake.wildcards.b}']
output_params = ['c', 'CpN', 'd', 'G', 'p', 'P', 'S']
for param in output_params:
    predicted_out.append(format_option_for_defiant_output(param))
predicted_out = '_'.join(predicted_out) + '.tsv'
move_command = f'mv {predicted_out} {out}'

# Execute shell command.
shell(
    "("
    "defiant "
    "-L {snakemake.wildcards.a},{snakemake.wildcards.b} "
    "-i {a_files} {b_files} "
    "{extra} "
    "{user_parameters} "
    "-cpu {snakemake.threads} && "
    "{move_command}"
    ")"
    "{log}"
)
