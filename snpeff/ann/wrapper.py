__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import sys
from os import path

from snakemake.shell import shell

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + snakemake.input[parameter]
    except AttributeError:
        return ''

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + snakemake.params[parameter]
    except AttributeError:
        return ''

# Define exception classes.
class RuleParameterException(Exception):
    pass

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
java_options = snakemake.params.get('java_options', '-Xmx4g')

genome_version = snakemake.params.get('genome_version', '')
if genome_version == '':
    raise RuleParameterException('You should specify a genome version that you used when creating your VCF file.')

no_statistics = snakemake.params.get('no_statistics', False)
no_statistics_option = ''
if no_statistics:
    no_statistics_option = '-noStats'

# Extract required arguments.
variants = snakemake.input[0]
if variants.endswith('vcf'):
    input_format = 'vcf'
elif variants.endswith('vcf.gz'):
    input_format = 'vcf'
    variants = variants[:-3]  # Trim tailing '.gz' from input.
elif vcf.endswith('bed'):
    input_format = 'bed'

output = snakemake.output[0]
if output.endswith('vcf'):
    output_format = 'vcf'
    pipe_command = '> %s' % output
elif output.endswith('vcf.gz'):
    output_format = 'vcf'
    pipe_command = '| gzip > %s' % output
elif output.endswith('bed'):
    output_format = 'bed'
    pipe_command = '> %s' % output


# Execute shell command.
shell(
    "("
    "snpEff "
    "{extra} "
    "{java_options} "
    "{no_statistics_option} "
    "-i {input_format} "
    "-o {output_format} "
    "{genome_version} "
    "{variants} "
    "{pipe_command}"
    ")"
    "{log}"
)
