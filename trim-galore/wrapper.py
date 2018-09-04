__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define exception classes.
class RuleParameterException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')
# `--fastqc` flag should not be included.
if '--fastqc' in extra:
    raise RuleParameterException('Please do not run trim_galore with --fastqc option.')

# Extract required arguments.
# Input should be single-ended or paired-ended.
# Single-end case.
if len(snakemake.input) == 1:
    read_command = snakemake.input[0]
# Paired-end case.
else:
    a, b = snakemake.input[0], snakemake.input[1]
    read_command = '--paired %s %s' % (a, b)

output_directory = path.dirname(snakemake.output[0])
output_directory_option = ''
if output_directory != '':
    output_directory_option = '-o ' + output_directory

# NOTE: For paired-end case, we rename output file *.read1_val_1.fq.gz into *.read1.trimmed.fastq.gz,
# and *.read2_val_2.fq.gz into *.read2.trimmed.fastq.gz
if len(snakemake.input) == 2:
    # Extract sample name.
    for output in snakemake.output:
        if output.endswith('.read1.trimmed.fastq.gz'):
            sample_name = path.basename(output)[:-23]

    raw_read1_file = path.join(output_directory, '%s.read1_val_1.fq.gz' % sample_name)
    renamed_read1_file = path.join(output_directory, '%s.read1.trimmed.fastq.gz' % sample_name)

    raw_read2_file = path.join(output_directory, '%s.read2_val_2.fq.gz' % sample_name)
    renamed_read2_file = path.join(output_directory, '%s.read2.trimmed.fastq.gz' % sample_name)

    rename_command = '&& mv %s %s ' \
                     '&& mv %s %s ' %\
                     (raw_read1_file, renamed_read1_file,
                      raw_read2_file, renamed_read2_file)

# NOTE: For single-end case, we rename *_trimmed.fastq.gz into *.trimmed.fastq.gz.
else:
    # Extract sample name from *.fastq.gz
    sample_name = path.basename(snakemake.input[0])[:-9]
    raw_read_file = path.join(output_directory, '%s_trimmed.fq.gz' % sample_name)
    renamed_read_file = path.join(output_directory, '%s.trimmed.fastq.gz' % sample_name)
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
    "{output_directory_option} "
    "{extra} "
    "{rename_command}"
    ")"
    "{log} "
)
