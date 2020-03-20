__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
bams = snakemake.input
merged_bam = snakemake.output

# Execute shell command.
shell(
    "("
    "sambamba merge "
    "{extra} "
    "-t {snakemake.threads} "
    "{merged_bam} "
    "{bams} "
    ") "
    "{log}"
)
