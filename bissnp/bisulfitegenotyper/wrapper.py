__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2021, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

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
user_parameters.append(optionify_params('stand_call_conf', '-stand_call_conf'))
user_parameters.append(optionify_params('min_mapping_quality_score', '-mmq'))
user_parameters.append(optionify_params('min_base_quality_score', '-mbq'))
user_parameters.append(optionify_params('output_modes', '--output_modes'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
bam = snakemake.input.bam
dbsnp = snakemake.input.dbsnp
reference = snakemake.input.reference

# Extract required outputs.
vcf_file_name1 = snakemake.output.vcf_file_name1
vcf_file_name2 = snakemake.output.vcf_file_name2

threads = snakemake.threads

# Execute shell command.
shell(
    "("
    "bis-snp "
    "-R {reference} "
    "-T BisulfiteGenotyper "
    "-I {bam} "
    "--dbsnp {dbsnp} "
    "-vfn1 {vcf_file_name1} "
    "-vfn2 {vcf_file_name2} "
    "{extra} "
    "{user_parameters} "
    "-nt {threads} "
    ")"
    "{log}"
)
