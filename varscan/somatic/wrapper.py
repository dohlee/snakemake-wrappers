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
snv_output = snakemake.output.snv
indel_output = snakemake.output.indel

output_prefix = snv_output[:snv_output.find('.snvs.varscan')]

# Extract optional parameters.
user_parameters = []
user_parameters.append(optionify_params('min_var_freq', '--min-var-freq'))
user_parameters.append(optionify_params('strand_filter', '--strand-filter'))
user_parameters = ' '.join(user_parameters)

wrapper_parameters = []
if not is_defined_by_user('--output-snp'):
    wrapper_parameters.append('--output-snp %s' % snv_output)
if not is_defined_by_user('--output-indel'):
    wrapper_parameters.append('--output-indel %s' % indel_output)
if not is_defined_by_user('--output-vcf'):
    wrapper_parameters.append('--output-vcf 1')
wrapper_parameters = ' '.join(wrapper_parameters)

# Pileup commands.
q_cutoff = snakemake.params.get('pileup_quality_cutoff', 20)
region_option = optionify_params('region', '-r')
normal_pileup_command = 'samtools mpileup -q %d -f %s %s %s' % (q_cutoff, reference, region_option, normal_bam)
tumor_pileup_command = 'samtools mpileup -q %d -f %s %s %s' % (q_cutoff, reference, region_option, tumor_bam)

# Execute shell command.
shell(
    "("
    "varscan somatic "
    "<({normal_pileup_command}) "
    "<({tumor_pileup_command}) "
    "{output_prefix} "
    "{extra} "
    "{user_parameters} "
    "{wrapper_parameters} "
    ") "
    "{log}"
)
