__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


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
user_parameters.append(optionify_params('kmerLen', '--kmerLen'))
user_parameters.append(optionify_params('gencode', '--gencode'))
user_parameters.append(optionify_params('keepDuplicates', '--keepDuplicates'))
user_parameters.append(optionify_params('filterSize', '--filterSize'))
user_parameters.append(optionify_params('tmpdir', '--tmpdir'))
user_parameters.append(optionify_params('sparse', '--sparse'))
user_parameters.append(optionify_params('decoys', '--decoys'))
user_parameters.append(optionify_params('type_', '--type'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required arguments.
reference = snakemake.input[0]
index = snakemake.output[0]

# Execute shell command.
shell(
    "("
    "salmon index "
    "--transcripts {reference} "
    "--index {index} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{user_parameters} "
    ")"
    "{log}"
)
