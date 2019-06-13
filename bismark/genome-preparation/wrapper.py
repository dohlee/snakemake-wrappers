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
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('verbose', '--verbose'))
user_parameters.append(optionify_params('path_to_aligner', '--path_to_aligner'))
user_parameters.append(optionify_params('bowtie2', '--bowtie2'))
user_parameters.append(optionify_params('hisat2', '--hisat2'))
user_parameters.append(optionify_params('single_fasta', '--single_fasta'))
user_parameters.append(optionify_params('genomic_composition', '--genomic_composition'))
user_parameters.append(optionify_params('slam', '--slam'))
user_parameters.append(optionify_params('large_index', '--large-index'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
reference_dir = snakemake.input

# Execute shell command.
shell(
    "("
    "bismark_genome_preparation "
    "{extra} "
    "{user_parameters} "
    "--parallel {snakemake.threads} "
    "{reference_dir} "
    ")"
    "{log}"
)
