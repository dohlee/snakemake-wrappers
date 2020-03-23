__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2020, Dohoon Lee"
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
if len(bams) > 1:
    # If more than one BAM files are given, merge them.
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
else:
    # Or if only one BAM file is given, just copy it.
    shell(
        "cp {bams} {merged_bam}"
    )
