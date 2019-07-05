__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import sys
import time
import uuid
from os import path

from snakemake.shell import shell

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
user_parameters.append(optionify_params('outputCallableRegions', '--outputCallabelRegions'))
user_parameters.append(optionify_params('indelCandidates', '--indelCandidates'))
user_parameters.append(optionify_params('forcedGT', '--forcedGT'))
user_parameters.append(optionify_params('exome', '--exome'))
user_parameters.append(optionify_params('callRegions', '--callRegions'))
user_parameters.append(optionify_params('noiseVcf', '--noiseVcf'))
user_parameters.append(optionify_params('scanSizeMb', '--scanSizeMb'))
user_parameters.append(optionify_params('region', '--region'))
user_parameters.append(optionify_params('callMemMb', '--callMemMb'))
user_parameters.append(optionify_params('retainTempFiles', '--retainTempFiles'))
user_parameters.append(optionify_params('disableEVS', '--disableEVS'))
user_parameters.append(optionify_params('reportEVSFeatures', '--reportEVSFeatures'))
user_parameters.append(optionify_params('snvScoringModelFile', '--snvScoringModelFile'))
user_parameters.append(optionify_params('indelScoringModelFile', '--indelScoringModelFile'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Randomly generate a name of temporary working directory.
run_directory = 'strelka_tmp_%s' % uuid.uuid4()

# Extract required input arguments.
reference = snakemake.input.reference
tumor_bam = snakemake.input.tumor
normal_bam = snakemake.input.normal

# Extract output
snv_output = snakemake.output.snv
indel_output = snakemake.output.indel

# Strelka is run in two Steps: configuration and execution.
configure_command = 'configureStrelkaSomaticWorkflow.py --tumorBam %s --normalBam %s --ref %s --runDir %s %s %s' % \
                    (tumor_bam, normal_bam, reference, run_directory, extra, user_parameters)
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
    "sleep 3 && "
    "{execution_command} && "
    "{move_command} && "
    "{clean_command}"
    ")"
    "{log}"
)
