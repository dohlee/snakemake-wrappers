__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
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
user_parameters.append(optionify_params('exome', '--exome'))
user_parameters.append(optionify_params('rna', '--rna'))
user_parameters.append(optionify_params('unstrandedRNA', '--unstrandedRNA'))
user_parameters.append(optionify_params('callRegions', '--callRegions'))
user_parameters.append(optionify_params('existingAlignStatsFile', '--existingAlignStatsFile'))
user_parameters.append(optionify_params('useExistingChromDepths', '--useExistingChromDepths'))
user_parameters.append(optionify_params('retainTempFiles', '--retainTempFiles'))
user_parameters.append(optionify_params('generateEvidenceBam', '--generateEvidenceBam'))
user_parameters.append(optionify_params('outputContig', '--outputContig'))
user_parameters.append(optionify_params('scanSizeMb', '--scanSizeMb'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Randomly generate a name of temporary working directory.
run_directory = 'manta_tmp_%s' % uuid.uuid4()

# Extract required input arguments.
reference = snakemake.input.reference
tumor_bam = snakemake.input.tumor
normal_bam = snakemake.input.normal

# Extract output
csi = snakemake.output.candidate_small_indels
csv = snakemake.output.candidate_sv
dsv = snakemake.output.diploid_sv
ssv = snakemake.output.somatic_sv
csi_idx = snakemake.output.candidate_small_indels_idx
csv_idx = snakemake.output.candidate_sv_idx
dsv_idx = snakemake.output.diploid_sv_idx
ssv_idx = snakemake.output.somatic_sv_idx
csi_orig = '%s/results/variants/candidateSmallIndels.vcf.gz' % run_directory
csv_orig = '%s/results/variants/candidateSV.vcf.gz' % run_directory
dsv_orig = '%s/results/variants/diploidSV.vcf.gz' % run_directory
ssv_orig = '%s/results/variants/somaticSV.vcf.gz' % run_directory

# Manta requires two steps: configuration and execution.
configure_command = 'configManta.py --tumorBam %s --normalBam %s --referenceFasta %s --runDir %s %s %s' % \
                    (tumor_bam, normal_bam, reference, run_directory, extra, user_parameters)
execution_command = '%s/runWorkflow.py -m local -j %d' % (run_directory, snakemake.threads)

# Move manta results to desired output directory.
move_command = ('mv %s %s ' * 8) % (csi_orig, csi, csv_orig, csv, dsv_orig, dsv, ssv_orig, ssv)
move_command = move_command.strip()

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
