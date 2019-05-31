__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import itertools
from os import path, makedirs

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define utility functions.
def get_prefix(f):
    return path.splitext(path.basename(f))[0]

def assert_directory_exists(directory):
    if not path.exists(directory):
        makedirs(directory)

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

# Define exception classes.
class RuleInputException(Exception):
    pass

class RuleParameterException(Exception):
    pass

class RuleOutputException(Exception):
    pass

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('color', '--color'))
user_parameters.append(optionify_params('noauto', '--noauto'))
user_parameters.append(optionify_params('packed', '--packed'))
user_parameters.append(optionify_params('bmax', '--bmax'))
user_parameters.append(optionify_params('bmaxdivn', '--bmaxdivn'))
user_parameters.append(optionify_params('dcv', '--dcv'))
user_parameters.append(optionify_params('nodc', '--nodc'))
user_parameters.append(optionify_params('noref', '--noref'))
user_parameters.append(optionify_params('justref', '--justref'))
user_parameters.append(optionify_params('offrate', '--offrate'))
user_parameters.append(optionify_params('ftabchars', '--ftabchars'))
user_parameters.append(optionify_params('ntoa', '--ntoa'))
user_parameters.append(optionify_params('seed', '--seed'))
user_parameters.append(optionify_params('quiet', '--quiet'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])


# Assert input and output have been correctly given.
if len(snakemake.input) != 1:
    raise RuleInputException('Please check your reference genome has been correctly given. It should be given as a single file.')

# Extract required inputs.
reference = snakemake.input.reference
prefix = get_prefix(reference)

# Extract required outputs.
index_dir = snakemake.output.index_dir
output = path.join(index_dir, prefix)

# Assert index directory exists.
assert_directory_exists(index_dir)

# Execute shell command.
shell(
    "("
    "bowtie-build "
    "{reference} "
    "{output} "
    "--threads {snakemake.threads}"
    "{extra} "
    "{user_parameters} "
    ") "
    "{log}"
)
