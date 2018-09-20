__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

def is_defined_by_user(*params):
    extra = snakemake.params.get('extra', '')
    for param in params:
        if param in extra:
            return True
    return False

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + str(snakemake.input[parameter])
    except AttributeError:
        return ''

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        return option + ' ' + str(snakemake.params[parameter])
    except AttributeError:
        return ''

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
tumor_bam = snakemake.input.tumor_bam
normal_bam = snakemake.input.normal_bam
reference = snakemake.input.reference

# Extract required outputs.
output = snakemake.output[0]

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params('min_alternate_fraction', '--min-alternate-fraction'))
user_parameters.append(optionify_params('min_alternate_count', '--min-alternate-count'))
user_parameters.append(optionify_params('cnv_map', '--cnv-map'))
user_parameters = ' '.join(user_parameters)

wrapper_parameters = []
if not is_defined_by_user('--pooled-discrete', '-J'):
    wrapper_parameters.append('-J')
if not is_defined_by_user('--pooled-continuous', '-K'):
    wrapper_parameters.append('-K')
if not is_defined_by_user('--report-genotype-likelihood-max'):
    wrapper_parameters.append('--report-genotype-likelihood-max')
if not is_defined_by_user('--allele-balance-priors-off'):
    wrapper_parameters.append('--allele-balance-priors-off')
if not is_defined_by_user('--min-repeat-entropy'):
    wrapper_parameters.append('--min-repeat-entropy 1')
if not is_defined_by_user('--no-partial-observations'):
    wrapper_parameters.append('--no-partial-observations')
wrapper_parameters = ' '.join(wrapper_parameters)

# Input bams option.
bams_option = '--bam %s --bam %s' % (tumor_bam, normal_bam)

# If user wants bcf file as an output, pipe the output through bcftools.
pipe_command = ''
if output.endswith('.bcf'):
    pipe_command += '| bcftools view -Ob -'

chunksize = snakemake.params.get('chunksize', 100000)

if snakemake.threads == 1:
    freebayes = 'freebayes'
else:
    freebayes = 'freebayes-parallel <(fasta_generate_regions.py %s.fai %d %d)' % (reference, chunksize, snakemake.threads)

# Execute shell command.
shell(
    "("
    "{freebayes} "
    "-f {reference} "
    "{extra} "
    "{user_parameters} "
    "{wrapper_parameters} "
    "{bams_option} "
    "{pipe_command} > "
    "{output} "
    ") "
    "{log}"
)
