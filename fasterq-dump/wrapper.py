__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define exception class.
class RuleOutputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
sra = snakemake.input

# Output should be one or two files
# *.read1.fastq.gz, *.read2.fastq.gz
if len(snakemake.output) not in [1, 2]:
    raise RuleOutputException('Please specify one(single-end) or two(paired-end) fastq files to the output.\n\n\
                               Hint: *.fastq.gz for single-end, *.read1.fastq.gz and *.read2.fastq.gz for paired-end.')
output_directory = path.dirname(snakemake.output[0]) or '.'

if len(snakemake.output) == 2:
    # Extract sample name.
    for output in snakemake.output:
        if output.endswith('.read1.fastq.gz'):
            sample_name = output[:-15]

    raw_read1_file = '%s.sra_1.fastq.gz' % sample_name
    renamed_read1_file = '%s.read1.fastq.gz' % sample_name

    raw_read2_file = '%s.sra_2.fastq.gz' % sample_name
    renamed_read2_file = '%s.read2.fastq.gz' % sample_name

    gzip_command = '&& pigz %s.sra_1.fastq --processes %d' \
                   '&& pigz %s.sra_2.fastq --processes %d' %\
                   (sample_name, snakemake.threads,
                    sample_name, snakemake.threads)

    rename_command = '&& mv %s %s ' \
                     '&& mv %s %s ' %\
                     (raw_read1_file, renamed_read1_file,
                      raw_read2_file, renamed_read2_file)

# NOTE: For single-end case, we rename *_pass.fastq.gz into *.fastq.gz.
else:
    # Extract sample name from *.fastq.gz
    sample_name = snakemake.output[0][:-9]
    raw_read_file = '%s.sra.fastq.gz' % sample_name
    renamed_read_file = '%s.fastq.gz' % sample_name

    gzip_command = '&& pigz %s.sra.fastq --processes %d' % (sample_name, snakemake.threads)
    rename_command = '&& mv %s %s' % (raw_read_file, renamed_read_file)

# NOTE: I fixed some recommended options for fastq-dump.
# Refer to: `fastq-dump-best-practice` in https://github.com/dohlee/bioinformatics-one-liners.
# Execute shell command.
shell(
    "("
    "fasterq-dump "
    "{sra} "
    "--threads {snakemake.threads} "
    "--outdir {output_directory} "
    "--split-3 --skip-technical "
    "{extra} "
    "{gzip_command} "
    "{rename_command} "
    ") "
    "{log} "
)
