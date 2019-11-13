__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
sorted_bam = snakemake.input[0]
flagstat = snakemake.output[0]

# Execute shell command.
shell(
    "("
    "sambamba flagstat "
    "{sorted_bam} "
    "{extra} "
    "-t {snakemake.threads} > "
    "{flagstat}"
    ") "
    "{log}"
)
