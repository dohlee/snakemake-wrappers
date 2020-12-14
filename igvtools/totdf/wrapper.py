__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2020, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Extract required inputs.
input = snakemake.input[0]

# Extract required outputs.
output = snakemake.output[0]

# Extract genome version parameter.
genome_version = snakemake.params['genome_version']

# Execute shell command.
shell(
    "("
    "igvtools toTDF "
    "{input} "
    "{output} "
    "{genome_version} "
    ")"
    "{log}"
)
