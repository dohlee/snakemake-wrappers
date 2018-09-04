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
read = snakemake.input.read  # Input read file in single-end(string) or paired-end(list).
reference = snakemake.input.reference  # Reference FASTA file.
output = snakemake.output  # Output bam/sam file.

# Command for single-end or paired-end reads.
read_command = '-a %s' % read if not isinstance(read, list) else '-a %s -b %s' % (read[0], read[1])

# Execute shell command.
shell(
    "("
    "bsmap "
    "{read_command} "
    "-d {reference} "
    "-o {output} "
    "-p {snakemake.threads} "
    "{extra} "
    ")"
    "{log}"
)
