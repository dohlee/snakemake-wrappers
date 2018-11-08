__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.input[parameter])
        if param:
            return option + ' ' + str(snakemake.input[parameter])
        else:
            return ''
    except AttributeError:
        return ''


def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.params[parameter])
        if param:
            return option + ' ' + str(snakemake.params[parameter])
        else:
            return ''
    except AttributeError:
        return ''


# Extract required inputs.
reference = snakemake.input.reference
assert not reference.endswith('.gz'), 'To use STAR, reference genome should be uncompressed.'

# Extract required outputs.
index_directory = snakemake.output.index_directory

mkdir_command = ''
if not os.path.exists(index_directory):
    mkdir_command = 'mkdir %s &&' % index_directory

# Extract parameters.
# Extract optional parameters.
extra = snakemake.params.get('extra', '')
sjdb_gtf_file = optionify_params('sjdb_gtf_file', '--sjdbGTFfile')
sjdb_overhang = optionify_params('sjdb_overhang', '--sjdbOverhang')
sjdb_gtf_chr_prefix = optionify_params('sjdb_gtf_chr_prefix', '--sjdbGTFchrPrefix')
sjdb_gtf_tag_exon_parent_transcript = optionify_params('sjdb_gtf_tag_exon_parent_transcript', '--sjdbGTFtagExonParentTranscript')

# Execute shell command.
shell(
    "("
    "{mkdir_command} "
    "STAR "
    "--runMode genomeGenerate "
    "--runThreadN {snakemake.threads} "
    "--genomeFastaFiles {reference} "
    "{extra} "
    "{sjdb_gtf_file} "
    "{sjdb_overhang} "
    "{sjdb_gtf_chr_prefix} "
    "{sjdb_gtf_tag_exon_parent_transcript} "
    "--genomeDir {index_directory}) "
    "{log}"
)
