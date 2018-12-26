__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import os

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

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
bams = snakemake.input.bams
targets = snakemake.input.targets
access = snakemake.input.access

# Extract required outputs.
output_target = snakemake.output.output_target
output_antitarget = snakemake.output.output_antitarget

# Generated target/antitarget bed file defaults to be in the working directory,
# so we should move them to the desired directory explicitly.
target_basename = os.path.splitext(os.path.basename(targets))[0]
generated_output_target = target_basename + '.target.bed'
generated_output_antitarget = target_basename + '.antitarget.bed'
move_command = f'mv {generated_output_target} {output_target} && mv {generated_output_antitarget} {output_antitarget}'

# Extract optional parameters.
user_parameters = []
user_parameters = ' '.join(user_parameters)

# Execute shell command.
shell(
    "("
    "cnvkit.py autobin {bams} "
    "--targets {targets} "
    "--access {access} "
    "{extra} && "
    "{move_command} "
    ") "
    "{log}"
)
