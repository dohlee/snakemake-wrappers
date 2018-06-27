__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
fastq = snakemake.input.fastq
fastq = [fastq] if isinstance(fastq, str) else fastq
assert len(fastq) <= 2, "Your sequencing read should be single- or paired-ended."
fastq = ' '.join(fastq)

index = snakemake.input.index
threads = snakemake.threads

output_directory = path.dirname(snakemake.output[0])

# Execute shell command.
shell(
    "("
    "kallisto quant "
    "-i {index} "
    "-o {output_directory} "
    "-t {threads} "
    "{extra} "
    "{fastq}"
    ")"
    "{log}"
)
