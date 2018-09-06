__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define exception classes.
class RuleInputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
fastq = snakemake.input.fastq
fastq = [fastq] if isinstance(fastq, str) else fastq
if len(fastq) > 2:
    raise RuleInputException('Your sequencing read should be single-read or paired-end.')
single_flag = '' if len(fastq) == 2 else '--single'
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
    "{single_flag} "
    "{extra} "
    "{fastq}"
    ")"
    "{log}"
)
