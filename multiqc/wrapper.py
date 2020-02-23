__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2020, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir, dirname
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
user_parameters.append(optionify_params('force', '--force'))
user_parameters.append(optionify_params('dirs', '--dirs'))
user_parameters.append(optionify_params('dirs_depth', '--dirs-depth'))
user_parameters.append(optionify_params('fullnames', '--fullnames'))
user_parameters.append(optionify_params('title', '--title'))
user_parameters.append(optionify_params('comment', '--comment'))
user_parameters.append(optionify_params('outdir', '--outdir'))
user_parameters.append(optionify_params('template', '--template'))
user_parameters.append(optionify_params('tag', '--tag'))
user_parameters.append(optionify_params('ignore', '--ignore'))
user_parameters.append(optionify_params('ignore_symlinks', '--ignore_symlinks'))
user_parameters.append(optionify_params('sample_names', '--sample-names'))
user_parameters.append(optionify_params('file_list', '--file-list'))
user_parameters.append(optionify_params('exclude', '--exclude'))
user_parameters.append(optionify_params('module', '--module'))
user_parameters.append(optionify_params('data_dir', '--data-dir'))
user_parameters.append(optionify_params('no_data_dir', '--no-data-dir'))
user_parameters.append(optionify_params('data_format', '--data-format'))
user_parameters.append(optionify_params('zip_data_dir', '--zip-data-dir'))
user_parameters.append(optionify_params('export', '--export'))
user_parameters.append(optionify_params('flat', '--flat'))
user_parameters.append(optionify_params('interactive', '--interactive'))
user_parameters.append(optionify_params('lint', '--lint'))
user_parameters.append(optionify_params('pdf', '--pdf'))
user_parameters.append(optionify_params('no_megaqc_upload', '--no-megaqc-upload'))
user_parameters.append(optionify_params('config', '--config'))
user_parameters.append(optionify_params('cl_config', '--cl-config'))
user_parameters.append(optionify_params('verbose', '--verbose'))
user_parameters.append(optionify_params('quiet', '--quiet'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
directory_name = dirname(snakemake.input[0])

# Extract required outputs.
output = snakemake.output[0]
output_directory = dirname(output)

# Always print out to stdout, and redirect it.
user_parameters.append(f'--outdir {output_directory}')
user_parameters.append('--filename stdout')

# Execute shell command.
shell(
    "("
    "multiqc {user_parameters} {directory_name} "
    "> {output}"
    ")"
    "{log}"
)
