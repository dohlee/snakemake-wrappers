__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define exception classes.
class RuleInputException(Exception):
    pass

class RuleParameterException(Exception):
    pass

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.input[parameter])
        if param:
            return option + ' ' + str(snakemake.input[parameter])
        else:
            return ''
    except AttributeError:
        return ''


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

# Extract parameters.
extra = snakemake.params.get('extra', '')
fragment_length = optionify_params('fragment_length', '-l')
standard_deviation = optionify_params('standard_deviation', '-s')

# Extract required arguments.
fastq = snakemake.input.fastq
fastq = [fastq] if isinstance(fastq, str) else fastq
if len(fastq) > 2:
    raise RuleInputException('Your sequencing read should be single-read or paired-end.')

single_flag = '' if len(fastq) == 2 else '--single'
if single_flag and (fragment_length == '' or standard_deviation == ''):
    raise RuleParameterException('Please provide fragment length(-l) and standard deviation(-s) parameter for single-end reads.')
fastq = ' '.join(fastq)

index = snakemake.input.index
threads = snakemake.threads

output_directory = path.dirname(snakemake.output[0])

# Execute shell command.
shell(
    "("
    "kallisto quant "
    "-i {index} "
    "-o {output_directory} "
    "-t {threads} "
    "{fragment_length} "
    "{standard_deviation} "
    "{single_flag} "
    "{extra} "
    "{fastq}"
    ")"
    "{log}"
)
