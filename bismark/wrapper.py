__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
from snakemake.shell import shell

# Define utility function.
def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        if str(snakemake.params[parameter]) == '':
            return ''
        if type(snakemake.params[parameter]) == bool:
            if snakemake.params[parameter]:
                return option
            else:
                return ''
        else:
            return option + ' ' + str(snakemake.params[parameter])
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('single_end', '--single_end'))
user_parameters.append(optionify_params('fastq', '--fastq'))
user_parameters.append(optionify_params('fasta', '--fasta'))
user_parameters.append(optionify_params('skip', '--skip'))
user_parameters.append(optionify_params('upto', '--upto'))
user_parameters.append(optionify_params('phred33_quals', '--phred33-quals'))
user_parameters.append(optionify_params('phred64_quals', '--phred64-quals'))
user_parameters.append(optionify_params('path_to_bowtie2', '--path_to_bowtie2'))
user_parameters.append(optionify_params('path_to_hisat2', '--path_to_hisat2'))
user_parameters.append(optionify_params('N', '-N'))
user_parameters.append(optionify_params('L', '-L'))
user_parameters.append(optionify_params('ignore_quals', '--ignore-quals'))
user_parameters.append(optionify_params('minins', '--minins'))
user_parameters.append(optionify_params('maxins', '--maxins'))
user_parameters.append(optionify_params('local', '--local'))
user_parameters.append(optionify_params('non_directional', '--non_directional'))
user_parameters.append(optionify_params('pbat', '--pbat'))
user_parameters.append(optionify_params('sam_no_hd', '--sam-no-hd'))
user_parameters.append(optionify_params('rg_tag', '--rg_tag'))
user_parameters.append(optionify_params('rg_id', '--rg_id'))
user_parameters.append(optionify_params('rg_sample', '--rg_sample'))
user_parameters.append(optionify_params('unmapped', '--unmapped'))
user_parameters.append(optionify_params('ambiguous', '--ambiguous'))
user_parameters.append(optionify_params('temp_dir', '--temp_dir'))
user_parameters.append(optionify_params('non_bs_mm', '--non_bs_mm'))
user_parameters.append(optionify_params('gzip', '--gzip'))
user_parameters.append(optionify_params('sam', '--sam'))
user_parameters.append(optionify_params('cram', '--cram'))
user_parameters.append(optionify_params('cram_ref', '--cram_ref'))
user_parameters.append(optionify_params('samtools_path', '--samtools_path'))
user_parameters.append(optionify_params('prefix', '--prefix'))
user_parameters.append(optionify_params('basename', '--basename'))
user_parameters.append(optionify_params('old_flag', '--old_flag'))
user_parameters.append(optionify_params('ambig_bam', '--ambig_bam'))
user_parameters.append(optionify_params('nucleotide_coverage', '--nucleotide_coverage'))
user_parameters.append(optionify_params('icpc', '--icpc'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
fastq = snakemake.input.fastq
reference_dir = snakemake.input.reference_dir

if len(fastq) == 2:
    input_params = '-1 %s -2 %s' % (fastq[0], fastq[1])
elif len(fastq) == 1:
    input_params = fastq[0]
else:
    raise ValueError('Please give 1 (single-read) or 2 (paired-end) fastq files for inputs of bismark.')

# Extract required outputs.
output_directory = path.dirname(snakemake.output[0])

# Scale threads since the actual number of threads in use will be about tripled (if using bowtie2).
threads = max(1, snakemake.threads // 3)

# Rename bismark outputs into e.g.
# 'result/se/{sample}/{sample}.bismark.bam'
# 'result/se/{sample}/{sample}.bismark_report.txt'
# or
# 'result/pe/{sample}/{sample}.bismark.bam'
# 'result/pe/{sample}/{sample}.bismark_report.txt'
basename = path.basename(fastq[0])[:-9] if fastq[0].endswith('.gz') else path.basename(fastq[0])[:-6]
if len(fastq) == 2:
    # Paired-end case.
    rename_command = '&& mv %s %s && mv %s %s' % (
        path.join(output_directory, basename + '_bismark_bt2_pe.bam'),
        snakemake.output[0],
        path.join(output_directory, basename + '_bismark_bt2_PE_report.txt'),
        snakemake.output[1],
    )
else:
    # Single-end case.
    rename_command = '&& mv %s %s && mv %s %s' % (
        path.join(output_directory, basename + '_bismark_bt2.bam'),
        snakemake.output[0],
        path.join(output_directory, basename + '_bismark_bt2_SE_report.txt'),
        snakemake.output[1],
    )

# Execute shell command.
shell(
    "("
    "bismark "
    "{extra} "
    "-o {output_directory} "
    "--parallel {threads} "
    "{user_parameters} "
    "{reference_dir} "
    "{input_params} "
    "{rename_command} "
    ")"
    "{log}"
)
