__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

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

# Assert input and output have been correctly given.
if len(snakemake.input) != 1:
    raise RuleInputException('Please check your reference genome has been correctly given. It should be given as a single file.')

if len(snakemake.output) != 5:
    raise RuleOutputException('bwa-mem generates 5 outputs, *.amb, *.ann, *.bwt, *.pac, and *.sa. Please check your specified output.')

# Extract required inputs.
reference = snakemake.input[0]
prefix = path.splitext(reference)[0]

algorithm = snakemake.params.get('algorithm', 'bwtsw')
# Assert the algorithm is 'is' or 'bwtsw'.
if algorithm not in ['is', 'bwtsw']:
    raise RuleParameterException('Algorithm should be "is" or "bwtsw"')

# Execute shell command.
shell(
    "("
    "bwa index "
    "-p {prefix} "
    "-a {algorithm} "
    "{extra} "
    "{reference} "
    ") "
    "{log}"
)
