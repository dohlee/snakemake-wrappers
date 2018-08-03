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
sra = snakemake.input

# Output should be one or three files
# *.read1.fastq.gz, *.read2.fastq.gz, *.fastq.gz
assert len(snakemake.output) in [1, 3], \
    'Please specify one(single-end) or three(paired-end) fastq.gz files to the output.\n\n\
    Hint: *.fastq.gz for single-end, *.read1.fastq.gz, *.read2.fastq.gz, *.fastq.gz for paired-end.'
output_directory = path.dirname(snakemake.output[0])

# If user wants output files to be gzipped, but did not specified --gzip option,
# kindly add --gzip option to extra options.
if all(f.endswith('.gz') for f in snakemake.output) and ('--gzip' not in extra):
    extra += ' --gzip'

# NOTE: I fixed some recommended options for fastq-dump.
# Refer to: `fastq-dump-best-practice` in https://github.com/dohlee/bioinformatics-one-liners.
# Execute shell command.
shell(
    "("
    "parallel-fastq-dump "
    "-s {sra} "
    "-t {snakemake.threads} "
    "-O {snakemake.output} "
    "--split-3 --skip-technical --readids --clip --read-filter pass "
    "{extra} "
    ")"
    "{log}"
)
