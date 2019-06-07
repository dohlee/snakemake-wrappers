__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
output = snakemake.output[0]

# Execute shell command.
shell(
    "("
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/bed/test.bed -qO {output}"
    ")"
)
