__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
from snakemake.shell import shell

# Define utility function.
def get_common_prefixes(strings):
    all_same = lambda x: all(x[0] == y for y in x)

    prefix_tuples = itertools.takewhile(all_same, zip(*strings))
    return ''.join(x[0] for x in prefix_tuples)

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
reads = snakemake.input.reads
if len(reads) == 2:
    read_command = '-1 %s -2 %s' % (reads[0], reads[1])
else:
    read_command = '%s' % (reads[0])

index_dir = snakemake.input.index_dir
prefix = get_common_prefixes(listdir(index_dir))
index_command = path.join(index_dir, prefix)

# Extract required outputs.
output = snakemake.output[0]
if output.endswith('.sorted.bam'):
    if 'mapq_cutoff' in snakemake.params:
        postprocess_command = '| samtools view -bS -q %d - | samtools sort > %s' % (snakemake.params.mapq_cutoff, output)
    else:
        postprocess_command = '| samtools view -bS - | samtools sort > %s' % (output)
else:
    if 'mapq_cutoff' in snakemake.params:
        postprocess_command = '| samtools view -bS -q %d - > %s' % (snakemake.params.mapq_cutoff, output)
    else:
        postprocess_command = '| samtools view -bS - > %s' % (output)

# Execute shell command.
shell(
    "("
    "bowtie "
    "{index_command} "
    "{read_command} "
    "{extra} "
    "--threads {snakemake.threads} "
    "{postprocess_command}) "
    "{log}"
)
