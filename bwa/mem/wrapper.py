__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


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

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('k', '-k'))
user_parameters.append(optionify_params('w', '-w'))
user_parameters.append(optionify_params('d', '-d'))
user_parameters.append(optionify_params('r', '-r'))
user_parameters.append(optionify_params('y', '-y'))
user_parameters.append(optionify_params('c', '-c'))
user_parameters.append(optionify_params('D', '-D'))
user_parameters.append(optionify_params('W', '-W'))
user_parameters.append(optionify_params('m', '-m'))
user_parameters.append(optionify_params('S', '-S'))
user_parameters.append(optionify_params('P', '-P'))
user_parameters.append(optionify_params('A', '-A'))
user_parameters.append(optionify_params('B', '-B'))
user_parameters.append(optionify_params('O', '-O'))
user_parameters.append(optionify_params('E', '-E'))
user_parameters.append(optionify_params('L', '-L'))
user_parameters.append(optionify_params('U', '-U'))
user_parameters.append(optionify_params('x', '-x'))
user_parameters.append(optionify_params('R', '-R'))
user_parameters.append(optionify_params('H', '-H'))
user_parameters.append(optionify_params('j', '-j'))
user_parameters.append(optionify_params('_5', '-5'))
user_parameters.append(optionify_params('q', '-q'))
user_parameters.append(optionify_params('K', '-K'))
user_parameters.append(optionify_params('v', '-v'))
user_parameters.append(optionify_params('T', '-T'))
user_parameters.append(optionify_params('h', '-h'))
user_parameters.append(optionify_params('a', '-a'))
user_parameters.append(optionify_params('C', '-C'))
user_parameters.append(optionify_params('V', '-V'))
user_parameters.append(optionify_params('M', '-M'))
user_parameters.append(optionify_params('I', '-I'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
reads = snakemake.input.reads
reference = snakemake.input.reference
db_prefix = path.splitext(reference)[0]

# Extract required outputs.
output = snakemake.output[0]

pipe_command = ''
# bwa-mem output defaults to sam. Convert sam to bam with samtools.
if output.endswith('.sorted.bam'):
    pipe_command = '| samtools view -Sb - | samtools sort'
elif output.endswith('.bam'):
    pipe_command = '| samtools view -Sb -'

# Execute shell command.
shell(
    "("
    "bwa mem "
    "{extra} "
    "{user_parameters} "
    "-t {snakemake.threads} "
    "{db_prefix} "
    "{reads} "
    "{pipe_command} > "
    "{output}) "
    "{log}"
)
