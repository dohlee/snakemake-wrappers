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
user_parameters.append(optionify_params('no_qualities', '--no-qualities'))
user_parameters.append(optionify_params('strandedness', '--strandedness'))
user_parameters.append(optionify_params('alignments', '--alignments'))
user_parameters.append(optionify_params('bowtie2', '--bowtie2'))
user_parameters.append(optionify_params('star', '--star'))
user_parameters.append(optionify_params('append_names', '--append-names'))
user_parameters.append(optionify_params('seed', '--seed'))
user_parameters.append(optionify_params('single_cell_prior', '--single_cell_prior'))
user_parameters.append(optionify_params('calc_pme', '--calc-pme'))
user_parameters.append(optionify_params('calc_ci', '--calc-ci'))
user_parameters.append(optionify_params('quiet', '--quiet'))
user_parameters.append(optionify_params('sort_bam_by_read_name', '--sort-bam-by-read-name'))
user_parameters.append(optionify_params('no_bam_output', '--no-bam-output'))
user_parameters.append(optionify_params('sampling_for_bam', '--sampling-for-bam'))
user_parameters.append(optionify_params('output_genome_bam', '--output-genome-bam'))
user_parameters.append(optionify_params('sort_bam_by_coordinate', '--sort-bam-by-coordinate'))
user_parameters.append(optionify_params('sort_bam_memory_per_thread', '--sort-bam-memory-per-thread'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required arguments.
reads = snakemake.input.reads
reference = snakemake.input.reference[:-15]  # Strip trailing '.transcripts.fa'.
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
    "{user_parameters} "
    "{extra} "
    ")"
    "{log}"
)
