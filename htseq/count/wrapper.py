__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

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

# Extract required inputs.
alignment = snakemake.input.alignment
annotation = snakemake.input.annotation

# Extract required outputs.
output = snakemake.output[0]

# Extract parameters.
# Extract optional parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('order', '--order'))
user_parameters.append(optionify_params('max_reads_in_buffer', '--max-reads-in-buffer'))
user_parameters.append(optionify_params('stranded', '--stranded'))
user_parameters.append(optionify_params('minaqual', '--minaqual'))
user_parameters.append(optionify_params('type', '--type'))
user_parameters.append(optionify_params('idattr', '--idattr'))
user_parameters.append(optionify_params('additional_attr', '--additional-attr'))
user_parameters.append(optionify_params('mode', '--mode'))
user_parameters.append(optionify_params('nonunique', '--nonunique'))
user_parameters.append(optionify_params('secondary_alignments', '--secondary-alignments'))
user_parameters.append(optionify_params('samout', '--samout'))
user_parameters.append(optionify_params('quiet', '--quiet'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

file_format_option = ''
if alignment.endswith('.bam') and '-f bam' not in extra:
    file_format_option = '-f bam'

# Execute shell command.
shell(
    "("
    "htseq-count "
    "{extra} "
    "{user_parameters} "
    "{file_format_option} "
    "{alignment} "
    "{annotation}"
    ") > {output} "
    "{log}"
)
