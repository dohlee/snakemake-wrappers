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
if isinstance(reads, str):
    reads = [reads]
if isinstance(mates, str):
    mates = [mates]
assert len(reads) >= 1, 'You should give at least one reads and mates. Given: %d' % len(reads)
assert len(reads) == len(mates), 'The number of reads and their pairs should match. Given: %d reads, %d mates' % (len(reads), len(mates))

reads, mates = ','.join(reads), ','.join(mates)

reference = snakemake.input.reference

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
    "bwameth.py "
    "{extra} "
    "-t {snakemake.threads} "
    "--reference {reference} "
    "{reads} "
    "{mates} "
    "{pipe_command} > "
    "{output}) "
    "{log}"
)
