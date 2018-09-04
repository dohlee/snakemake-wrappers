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

class StdioException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')
# `--stdin` and `--stdout` jflag should not be included.
if '--stdin' in extra or '--stdout' in extra:
    raise StdioException('Please do not run fastp with --stdin or --stdout option for now.')

# Extract required arguments.
# Input should be single-ended or paired-ended.
if len(snakemake.input) not in [1, 2]:
    raise RuleInputException('Input should be single-read or paired-end.')
# Input and output length should be the same.
if len(snakemake.input) != len(snakemake.output):
    raise RuleInputException('The number of input(%d) and output(%d) read files should be the same.' % (len(snakemake.input), len(snakemake.output)))
# Single-end case.
if len(snakemake.input) == 1:
    read_command = '-i %s' % snakemake.input[0]
    out_command = '-o %s' % snakemake.output[0]
# Paired-end case.
else:
    a, b = snakemake.input[0], snakemake.input[1]
    read_command = '-i %s -I %s' % (a, b)

    a, b = snakemake.output[0], snakemake.output[1]
    out_command = '-o %s -O %s' % (a, b)

shell(
    "("
    "fastp "
    "{read_command} "
    "{out_command} "
    "-w {snakemake.threads} "
    "{extra} "
    ")"
    "{log} "
)
