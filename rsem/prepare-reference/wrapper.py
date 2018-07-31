__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
fasta = snakemake.input.fasta
gtf = snakemake.input.get('gtf', None)
gff3 = snakemake.input.get('gff3', None)

assert not ((gtf is not None) and (gff3 is not None)), 'You cannot provide both GTF and GFF3 files.'
annotation_option = ''
if gtf is not None:
    annotation_option = '--gtf {gtf}'
elif gff3 is not None:
    annotation_option = '--gff3 {gff3}'

output_prefix = snakemake.output[0].rstrip('.transcripts.fa')
threads = snakemake.threads

# Execute shell command.
shell(
    "("
    "rsem-prepare-reference "
    "{annotation_option} "
    "{extra} "
    "{fasta} "
    "{output_prefix} "
    ")"
    "{log}"
)
