__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
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

def is_bigwig(f):
    ext = path.splitext(f)[1]
    return ext in ['.bw', '.bigwig', '.bigWig']

def is_bedgraph(f):
    ext = path.splitext(f)[1]
    return ext in ['.bedgraph', '.bedGraph']

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
treatment_bam = snakemake.input.treatment_bam
control_bam = snakemake.input.control_bam

# Extract required outputs.
output = snakemake.output.output
# Determine output file format.
if is_bigwig(output):
    out_format = 'bigwig'
elif is_bedgraph(output):
    out_format = 'bedgraph'
else:
    raise TypeError('Output extension should be one of ["bw", "bigwig", "bigWig", "bedgraph", "bedGraph"]')

# Extract user parameters.
user_parameters = []
user_parameters.append(optionify_params('scale_factors_method', '--scaleFactorsMethod'))
user_parameters.append(optionify_params('operation', '--operation'))
user_parameters.append(optionify_params('bin_size', '--binSize'))
user_parameters.append(optionify_params('region', '--region'))
user_parameters.append(optionify_params('effective_genome_size', '--effectiveGenomeSize'))
user_parameters.append(optionify_params('normalize_using', '--normalizeUsing'))
user_parameters.append(optionify_params('skip_non_covered_regions', '--skipNonCoveredRegions'))
user_parameters.append(optionify_params('smooth_length', '--smoothLength'))
user_parameters.append(optionify_params('extend_reads', '--extendReads'))
user_parameters = ' '.join(user_parameters)
    
# Execute shell command.
shell(
    "("
    "bamCompare "
    "--bamfile1 {treatment_bam} "
    "--bamfile2 {control_bam} "
    "--outFileName {output} "
    "--outFileFormat {out_format} "
    "--numberOfProcessors {snakemake.threads} "
    "{user_parameters} "
    "{extra} "
    ") "
    "{log}"
)
