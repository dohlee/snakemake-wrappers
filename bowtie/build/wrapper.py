__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import itertools
from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define utility functions.
def get_prefix(f):
    return path.splitext(path.basename(f))[0]

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

# Extract required inputs.
reference = snakemake.input.reference
prefix = get_prefix(reference)

# Extract required outputs.
index_dir = snakemake.output.index_dir
output = path.join(index_dir, prefix)

# Execute shell command.
shell(
    "("
    "bowtie-build "
    "{reference} "
    "{output} "
    "--threads {snakemake.threads}"
    "{extra} "
    ") "
    "{log}"
)
