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

# Output should be one or two files
# *.read1.fastq.gz, *.read2.fastq.gz
assert len(snakemake.output) in [1, 2], \
    'Please specify one(single-end) or three(paired-end) fastq.gz files to the output.\n\n\
    Hint: *.fastq.gz for single-end, *.read1.fastq.gz, *.read2.fastq.gz for paired-end.'
output_directory = path.dirname(snakemake.output[0]) or '.'

if len(snakemake.output) == 2:
    # Extract sample name.
    for output in snakemake.output:
        if output.endswith('.read1.fastq.gz'):
            sample_name = output[:-15]

    raw_read1_file = '%s_pass_1.fastq.gz' % sample_name
    renamed_read1_file = '%s.read1.fastq.gz' % sample_name

    raw_read2_file = '%s_pass_2.fastq.gz' % sample_name
    renamed_read2_file = '%s.read2.fastq.gz' % sample_name

    rename_command = '&& mv %s %s ' \
                     '&& mv %s %s ' %\
                     (raw_read1_file, renamed_read1_file,
                      raw_read2_file, renamed_read2_file)

# NOTE: For single-end case, we rename *_pass.fastq.gz into *.fastq.gz.
else:
    # Extract sample name from *.fastq.gz
    sample_name = snakemake.output[0][:-9]
    raw_read_file = '%s_pass.fastq.gz' % sample_name
    renamed_read_file = '%s.fastq.gz' % sample_name
    rename_command = '&& mv %s %s' % (raw_read_file, renamed_read_file)

# If user wants output files to be gzipped, but did not specified --gzip option,
# kindly add --gzip option to extra options.
if all(f.endswith('.gz') for f in snakemake.output) and ('--gzip' not in extra):
    extra += '--gzip'

# NOTE: I fixed some recommended options for fastq-dump.
# Refer to: `fastq-dump-best-practice` in https://github.com/dohlee/bioinformatics-one-liners.
# Execute shell command.
shell(
    "("
    "parallel-fastq-dump "
    "-s {sra} "
    "-t {snakemake.threads} "
    "-O {output_directory} "
    "--split-3 --skip-technical --clip --read-filter pass "
    "{extra} "
    "{rename_command} "
    ") "
    "{log} "
)
