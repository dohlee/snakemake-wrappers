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
sorted_bam = snakemake.input[0]
indexed_bam = snakemake.output[0]

# Execute shell command.
shell(
    "("
    "sambamba index "
    "{extra} "
    "-t {snakemake.threads} "
    "{sorted_bam} "
    "{indexed_bam} "
    ") "
    "{log}"
)
