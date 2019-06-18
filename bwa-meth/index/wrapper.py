__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Execute shell command.
shell(
    "("
    "bwameth.py index "
    "{snakemake.input[0]} "
    ")"
    "{log}"
)
