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

# Extract logself.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
input = snakemake.input[0]

if input.endswith('.gff.gz'):
    preset_option = '-p gff'
elif input.endswith('.bed.gz'):
    preset_option = '-p bed'
elif input.endswith('.vcf.gz'):
    preset_option = '-p vcf'

# Execute shell command.
shell(
    "("
    "tabix "
    "{extra} "
    "{preset_option} "
    "{input} "
    ") "
    "{log}"
)
