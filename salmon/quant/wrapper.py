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
library_type = snakemake.params.get('library_type', 'A')

# Extract required arguments.
fastq = snakemake.input.fastq
fastq = [fastq] if isinstance(fastq, str) else fastq
assert len(fastq) <= 2, "Your sequencing read should be single- or paired-ended."

# Single-end command
if len(fastq) == 1:
    read_command = '-r ' + fastq[0]
# Paired-end command
else:
    read_command = '-1 ' + fastq[0] + ' -2 ' + fastq[1]

index = snakemake.input.index
threads = snakemake.threads

quant = snakemake.output.quant
lib = snakemake.output.lib
output_directory = path.dirname(snakemake.output.quant)

# Execute shell command.
shell(
    "("
    "salmon quant "
    "-i {index} "
    "{read_command} "
    "-l {library_type} "
    "-o {output_directory} "
    "-p {threads} "
    "{extra} "
    ")"
    "{log}"
)
