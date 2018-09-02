__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

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
if isinstance(reference, str):
    reference = [reference]
reference = ','.join(reference)
annotation = snakemake.input.annotation

# Extract required outputs.
output_prefix = snakemake.output[0][:-6]

# Extract parameters.
# Extract optional parameters.
extra = snakemake.params.get('extra', '')
snp = optionify_params('snp', '--mode')
haplotype = optionify_params('haplotype', '--stranded')
splice_site = optionify_params('splice_site', '--ss')
exon = optionify_params('exon', '--exon')

# Execute shell command.
shell(
    "("
    "hisat2-build "
    "-p {snakemake.threads} "
    "{extra} "
    "{snp} "
    "{haplotype} "
    "{splice_site} "
    "{exon} "
    "{reference} "
    "{annotation} "
    ") "
    "{log}"
)
