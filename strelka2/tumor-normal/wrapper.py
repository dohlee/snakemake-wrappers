__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import sys
import time
from os import path

from snakemake.shell import shell


def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + snakemake.params[parameter]
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
# Extract parameters.
extra = snakemake.params.get('extra', '')
extra += optionify_params('region', '--region')
extra += optionify_params('call_regions', '--callRegions')

# Generate a name of temporary working directory.
run_directory = 'strelka_tmp_%d' % int(time.time() * 100)

# Extract required input arguments.
reference = snakemake.input.reference
tumor_bam = snakemake.input.tumor
normal_bam = snakemake.input.normal

# Extract output
snv_output = snakemake.output.snv
indel_output = snakemake.output.indel

# Strelka is run in two Steps: configuration and execution.
configure_command = 'configureStrelkaSomaticWorkflow.py --tumorBam %s --normalBam %s --ref %s --runDir %s %s' % \
                    (tumor_bam, normal_bam, reference, run_directory, extra)
execution_command = '%s/runWorkflow.py -m local -j %d' % (run_directory, snakemake.threads)

# Move strelka results to desired output directory.
move_command = 'mv %s/results/variants/somatic.snvs.vcf.gz %s && mv %s/results/variants/somatic.indels.vcf.gz %s' % \
                (run_directory, snv_output, run_directory, indel_output)

# Clean temporary directories.
clean_command = 'rm -r ' + run_directory

# Execute shell command.
shell(
    "("
    "{configure_command} && "
    "{execution_command} && "
    "{move_command} && "
    "{clean_command}"
    ")"
    "{log}"
)
