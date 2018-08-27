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

# Extract required inputs.
reads = snakemake.input.reads
mates = snakemake.input.get('mates', '')
reference = snakemake.input.reference

db_prefix = path.splitext(reference)[0]

# If -p option is present, 'mates' will be ignored and bwa-mem will assume 2i and 2i+1-th read files are paired.
input_is_interleaved = '-p' in extra
if input_is_interleaved and mates != '':  # Kindly reformat input data for mistaken input.
    reads = [file for pair in zip(read, mate) for file in pair]  # Interleave two lists!
    mates = ''

# Extract required outputs.
output = snakemake.output[0]

pipe_command = ''
# bwa-mem output defaults to sam. Convert sam to bam with samtools.
# TODO: Use sambamba with appropriate number of threads.
if output.endswith('.bam'):
    pipe_command = '| samtools view -Sb -'

# Execute shell command.
shell(
    "("
    "bwa mem "
    "{extra} "
    "-t {snakemake.threads} "
    "{db_prefix} "
    "{reads} "
    "{mates} "
    ") "
    "{pipe_command} > "
    "{output} "
    "{log}"
)
