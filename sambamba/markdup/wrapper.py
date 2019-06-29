__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

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

# Extract required arguments.
sorted_bam = snakemake.input[0]
duplicates_marked_bam = snakemake.output[0]

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('remove_duplicates', '--remove-duplicates'))
user_parameters.append(optionify_params('compression_level', '--compression-level'))
user_parameters.append(optionify_params('show_progress', '--show-progress'))
user_parameters.append(optionify_params('hash_table_size', '--hash-table-size'))
user_parameters.append(optionify_params('overflow_list_size', '--overflow-list-size'))
user_parameters.append(optionify_params('sort_buffer_size', '--sort-buffer-size'))
user_parameters.append(optionify_params('io_buffer_size', '--io-buffer-size'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Execute shell command.
shell(
    "("
    "sambamba markdup "
    "{extra} "
    "{user_parameters} "
    "-t {snakemake.threads} "
    "{sorted_bam} "
    "{duplicates_marked_bam} "
    ") "
    "{log}"
)
