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
# `--fastqc` flag should not be included.
assert '--fastqc' not in extra, "Don't run trim_galore with --fastqc option."

# Extract required arguments.
# Input should be single-ended or paired-ended.
# Single-end case.
if len(snakemake.input) == 1:
    read_command = snakemake.input[0]
# Paired-end case.
else:
    read_command = '--paired {snakemake.input[0]} {snakemake.input[1]}'

output_directory = path.dirname(snakemake.output[0])
# NOTE: For paired-end case, we rename output file *.read1_val_1.fq.gz into *.read1.trimmed.fastq.gz,
# and *.read2_val_2.fq.gz into *.read2.trimmed.fastq.gz
if len(snakemake.input) == 2:
    # Extract sample name.
    for output in snakemake.output:
        if output.endswith('.read1.fastq.gz'):
            sample_name = output[:-15]

    raw_read1_file = '%s.read1_val_1.fq.gz' % sample_name
    renamed_read1_file = '%s.read1.trimmed.fastq.gz' % sample_name

    raw_read2_file = '%s.read2_val_2.fq.gz' % sample_name
    renamed_read2_file = '%s.read2.trimmed.fastq.gz' % sample_name

    rename_command = '&& mv %s %s ' \
                     '&& mv %s %s ' %\
                     (raw_read1_file, renamed_read1_file,
                      raw_read2_file, renamed_read2_file)

# NOTE: For single-end case, we rename *_trimmed.fastq.gz into *.trimmed.fastq.gz.
else:
    # Extract sample name from *.fastq.gz
    sample_name = snakemake.input[0][:-9]
    raw_read_file = '%s_trimmed.fq.gz' % sample_name
    renamed_read_file = '%s.trimmed.fastq.gz' % sample_name
    rename_command = '&& mv %s %s' % (raw_read_file, renamed_read_file)

# If user wants output files to be gzipped, but did not specified --gzip option,
# kindly add --gzip option to extra options.
if all(f.endswith('.gz') for f in snakemake.output) and ('--gzip' not in extra):
    extra += ' --gzip'

# NOTE: I fixed some recommended options for fastq-dump.
# Refer to: `fastq-dump-best-practice` in https://github.com/dohlee/bioinformatics-one-liners.
# Execute shell command.
shell(
    "("
    "trim_galore "
    "{extra} "
    "{read_command} "
    "-o {output_directory} "
    "{extra} "
    ")"
    "{log} "
    "{rename_command}"
)
