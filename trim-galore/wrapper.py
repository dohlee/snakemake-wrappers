__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from os import path
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
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('quality', '--quality'))
user_parameters.append(optionify_params('phred33', '--phred33'))
user_parameters.append(optionify_params('phred64', '--phred64'))
user_parameters.append(optionify_params('fastqc', '--fastqc'))
user_parameters.append(optionify_params('fastqc_args', '--fastqc_args'))
user_parameters.append(optionify_params('adapter', '--adapter'))
user_parameters.append(optionify_params('adapter2', '--adapter2'))
user_parameters.append(optionify_params('illumina', '--illumina'))
user_parameters.append(optionify_params('nextera', '--nextera'))
user_parameters.append(optionify_params('small_rna', '--small_rna'))
user_parameters.append(optionify_params('max_length', '--max_length'))
user_parameters.append(optionify_params('stringency', '--stringency'))
user_parameters.append(optionify_params('e', '-e'))
user_parameters.append(optionify_params('length', '--length'))
user_parameters.append(optionify_params('max_n', '--max_n'))
user_parameters.append(optionify_params('trim_n', '--trim-n'))
user_parameters.append(optionify_params('no_report_file', '--no_report_file'))
user_parameters.append(optionify_params('suppress_warn', '--suppress_warn'))
user_parameters.append(optionify_params('clip_R1', '--clip_R1'))
user_parameters.append(optionify_params('clip_R2', '--clip_R2'))
user_parameters.append(optionify_params('three_prime_clip_R1', '--three_prime_clip_R1'))
user_parameters.append(optionify_params('three_prime_clip_R2', '--three_prime_clip_R2'))
user_parameters.append(optionify_params('nextseq', '--nextseq'))
user_parameters.append(optionify_params('basename', '--basename'))
user_parameters.append(optionify_params('rrbs', '--rrbs'))
user_parameters.append(optionify_params('non_directional', '--non_directional'))
user_parameters.append(optionify_params('keep', '--keep'))
user_parameters.append(optionify_params('trim1', '--trim1'))
user_parameters.append(optionify_params('retain_unpaired', '--retain_unpaired'))
user_parameters.append(optionify_params('length_1', '--length_1'))
user_parameters.append(optionify_params('length_2', '--length_2'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
if len(snakemake.input) == 1:
    read_params = snakemake.input[0]
elif len(snakemake.input) == 2:
    read_params = '--paired %s %s' % (snakemake.input[0], snakemake.input[1])
else:
    raise ValueError('Trim Galore! wrapper gets one or two files as input. Given: %d' % len(snakemake.input))

# Extract required outputs.
output_dir = path.dirname(snakemake.output[0])
output_dir_params = ''
if output_dir_params != '':
    output_dir_params = '-o ' + output_dir

# For paired-end case, we rename output file *.read1_val_1.fq.gz into *.read1.trimmed.fastq.gz,
# and *.read2_val_2.fq.gz into *.read2.trimmed.fastq.gz
if len(snakemake.input) == 2:
    # Extract sample name.
    for output in snakemake.output:
        if output.endswith('.read1.trimmed.fastq.gz'):
            sample_name = path.basename(output)[:-23]

    raw_read1_file = path.join(output_dir, '%s.read1_val1_fq.gz' % sample_name)
    renamed_read1_file = path.join(output_dir, '%s.read1.trimmed.fastq.gz' % sample_name)
    raw_read2_file = path.join(output_dir, '%s.read2_val_2.fq.gz' % sample_name)
    renamed_read2_file = path.join(output_dir, '%s.read2.trimmed.fastq.gz' % sample_name)

    rename_command = '&& sleep 5 && mv %s %s && mv %s %s' % (raw_read1_file, renamed_read1_file, raw_read2_file, renamed_read2_file)

# For single-read case, we rename *_trimmed.fq.gz into *.trimmed.fastq.gz.
else:
    # Extract sample name from *.fastq.gz
    sample_name = path.basename(snakemake.input[0])[:-9]
    raw_read_file = path.join(output_dir, '%s_trimmed.fq.gz' % sample_name)
    renamed_read_file = path.join(output_dir, '%s.trimmed.fastq.gz' % sample_name)
    rename_command = '&& mv %s %s' % (raw_read_file, renamed_read_file)

# If user wants output files to be gzipped, add --gzip option to extra options.
if all(f.endswith('.gz') for f in snakemake.output) and ('--gzip' not in extra):
    extra += ' --gzip'

# Execute shell command.
shell(
    "("
    "trim_galore "
    "{read_params} "
    "{extra} "
    "{user_parameters} "
    "--cores {snakemake.threads} "
    "{rename_command} "
    ") "
    "{log}"
)
