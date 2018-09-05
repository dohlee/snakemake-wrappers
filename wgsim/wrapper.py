__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

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

# Define exception classes.
class RuleInputException(Exception):
    pass

class RuleOutputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
if len(snakemake.input) > 1:
    raise RuleInputException('Please provide single reference genome file.')
reference = snakemake.input[0]

# Extract required outputs.
if len(snakemake.output) > 2:
    raise RuleOutputException('Simulated output should be single-read or paired-end.')
output = snakemake.output

# Extract parameters.
random_seed = optionify_params('random_seed', '-S')
N = optionify_params('N', '-N')
error_rate = optionify_params('error_rate', '-e')
mutation_rate = optionify_params('mutation_rate', '-r')
read_distance = optionify_params('read_distance', '-d')
standard_deviation = optionify_params('standard_deviation', '-s')
read1_length = optionify_params('read1_length', '-1')
read2_length = optionify_params('read2_length', '-2')
indel_fraction = optionify_params('indel_fraction', '-R')
indel_extension_probability = optionify_params('indel_extension_probability', '-X')
discard_ambiguous_bases_fraction = optionify_params('discard_ambiguous_bases_fraction', '-A')
haplotype_mode = ''
if snakemake.params.get('haplotype_mode', False):
    haplotype_mode = '-h'

# Execute shell command.
shell(
    "("
    "wgsim "
    "{random_seed} "
    "{N} "
    "{error_rate} "
    "{mutation_rate} "
    "{read_distance} "
    "{standard_deviation} "
    "{read1_length} "
    "{read2_length} "
    "{indel_fraction} "
    "{indel_extension_probability} "
    "{discard_ambiguous_bases_fraction} "
    "{haplotype_mode} "
    "{extra} "
    "{reference} "
    "{output}) "
    "{log}"
)
