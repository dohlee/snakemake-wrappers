__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from snakemake.shell import shell


# Helper function for finding common prefix.
def all_same(x):
    return all(x[0] == y for y in x)


def get_prefix_of_strings(strings):
    char_tuples = zip(*strings)
    prefix_tuples = itertools.takewhile(all_same, char_tuples)
    return ''.join(x[0] for x in prefix_tuples).strip('.')


# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
reads = snakemake.input.reads
reference = snakemake.input.reference
output_prefix = get_prefix_of_strings(snakemake.output)
threads = snakemake.threads

# Execute shell command.
shell(
    "("
    "rsem-calculate-expression "
    "--paired-end {reads} "
    "-p {threads} "
    "{reference} "
    "{output_prefix} "
    "{extra} "
    ")"
    "{log}"
)
