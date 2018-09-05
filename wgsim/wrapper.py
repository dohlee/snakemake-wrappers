__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define exception classes.
class RuleInputException(Exception):
    pass

class RuleOutputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
if len(snakemake.input) > 1:
    raise RuleInputException('Please provide single reference genome file.')
reference = snakemake.input[0]

# Extract required outputs.
if len(snakemake.output) > 2:
    raise RuleOutputException('Simulated output should be single-read or paired-end.')
output = snakemake.output

# Execute shell command.
shell(
    "("
    "wgsim "
    "{extra} "
    "-t {snakemake.threads} "
    "--reference {reference} "
    "{reads} "
    "{mates} "
    "{pipe_command} > "
    "{output}) "
    "{log}"
)
