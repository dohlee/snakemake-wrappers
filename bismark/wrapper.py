__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import os.path as path
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Define exception classes.
class GenomePreparationNotCompletedException(Exception):
    pass

class RuleInputException(Exception):
    pass

class RuleOutputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
fastq = snakemake.input.fastq
fastq = [fastq] if isinstance(fastq, str) else fastq
if len(fastq) > 2:
    raise RuleInputException('Your sequencing read should be single-read or paired-ended.')
read_command = '-1 %s -2 %s' % (fastq[0], fastq[1]) if len(fastq) == 2 else fastq[0]

# Path to the directory containing the unmodified reference genome,
# as well as the subfolders created by the `bismark_genome_preparation` script.
# Refer to: https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
reference_dir = snakemake.input.reference_dir

# Ensure `bismark_genome_preparation` script had been run.
bisulfite_genome_dir = snakemake.input.bisulfite_genome_dir
if path.join(reference_dir, 'Bisulfite_Genome') != bisulfite_genome_dir:
    raise GenomePreparationNotCompletedException('Please check that bismark_genome_preparation has been successfully finised.')

# Determine the number of threads.
# Since a typical Bismark run with 1 thread already uses about 2 (with --bowtie1) threads,
# and 3 (with --bowtie2 which is default option), this wrapper divides the number of
# user-defined threads with 2 or 3, depending on the bowtie option.
if '--bowtie1' in extra:
    threads = max(1, snakemake.threads // 2)
else:
    threads = max(1, snakemake.threads // 3)

# NOTE:
# Single-end Bismark call (with bowtie2) will produce two output files:
# 1. {fastq-filename}_bismark_bt2.bam
# 2. {fastq-filename}_bismark_bt2_SE_report.txt
# Paired-end Bismark call (with bowtie2) will produce two output files:
# 1. {fastq-read1-filename}_bismark_bt2_pe.bam
# 2. {fastq-read1-filename}_bismark_bt2_PE_report.txt

# Check there are two output file specified.
if len(snakemake.output) != 2:
    if len(fastq) == 1:
        raise RuleOutputException('Bismark has two outputs (with bowtie2): *_bismark_bt2.bam and *_bismark_bt2_SE_report.txt')
    else:
        raise RuleOutputException('Bismark has two outputs (with bowtie2): *_bismark_bt2_pe.bam and *_bismark_bt2_PE_report.txt')

output_directory = path.dirname(snakemake.output[0])

# Rename bismark outputs into
# 'result/{sample}/{sample}.bismark.bam',
# 'result/{sample}/{sample}.bismark_report.txt'
basename = fastq[0][:-9] if fastq[0].endswith('.gz') else fastq[0][:-6]
if len(fastq) == 2:
    # Paired-end case.
    rename_command = '&& mv %s %s && mv %s %s' % (
        basename + '_bismark_bt2_pe.bam', snakemake.output[0], basename + '_bismark_bt2_PE_report.txt', snakemake.output[1]
    )
else:
    # Single-end case.
    rename_command = '&& mv %s %s && mv %s %s' % (
        basename + '_bismark_bt2.bam', snakemake.output[0], basename + '_bismark_bt2_SE_report.txt', snakemake.output[1]
    )

# Execute shell command.
shell(
    "("
    "bismark "
    "{extra} "
    "-o {output_directory} "
    "--multicore {threads} "
    "{reference_dir} "
    "{read_command} "
    "{rename_command} "
    ")"
    "{log}"
)
