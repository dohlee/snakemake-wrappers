__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

import itertools

from os import path, listdir
from snakemake.shell import shell

# Define utility function.
def get_common_prefixes(strings):
    all_same = lambda x: all(x[0] == y for y in x)

    prefix_tuples = itertools.takewhile(all_same, zip(*strings))
    return ''.join(x[0] for x in prefix_tuples).strip('.')

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

# Extract required inputs.
treatment = snakemake.input.treatment
input_command = '-t %s' % treatment

# Extract optional inputs.
control = snakemake.input.control
if control:
    input_command += ' -c %s' % control

# Extract required outputs.
output_dir = path.dirname(snakemake.output.peak)
output_command = '--outdir %s' % output_dir

# Extract common prefixes of output files.
common_prefix = get_common_prefixes([path.basename(f) for f in snakemake.output])[:-1]
name_command = '-n %s' % common_prefix

# Extract user parameters.
user_parameters = []
user_parameters.append(optionify_params('gsize', '--gsize'))
user_parameters.append(optionify_params('keep_dup', '--keep-dup'))
user_parameters.append(optionify_params('buffer_size', '--buffer-size'))
user_parameters.append(optionify_params('name', '--name'))
user_parameters.append(optionify_params('bdg', '--bdg'))
user_parameters.append(optionify_params('trackline', '--trackline'))
user_parameters.append(optionify_params('SPMR', '--SPMR'))
user_parameters.append(optionify_params('tsize', '--tsize'))
user_parameters.append(optionify_params('bw', '--bw'))
user_parameters.append(optionify_params('mfold', '--mfold'))
user_parameters.append(optionify_params('fix_bimodal', '--fix-bimodal'))
user_parameters.append(optionify_params('nomodel', '--nomodel'))
user_parameters.append(optionify_params('shift', '--shift'))
user_parameters.append(optionify_params('extsize', '--extsize'))
user_parameters.append(optionify_params('qvalue', '--qvalue'))
user_parameters.append(optionify_params('pvalue', '--pvalue'))
user_parameters.append(optionify_params('to_large', '--to-large'))
user_parameters.append(optionify_params('ratio', '--ratio'))
user_parameters.append(optionify_params('down_sample', '--down-sample'))
user_parameters.append(optionify_params('seed', '--seed'))
user_parameters.append(optionify_params('tempdir', '--tempdir'))
user_parameters.append(optionify_params('nolambda', '--nolambda'))
user_parameters.append(optionify_params('slocal', '--slocal'))
user_parameters.append(optionify_params('llocal', '--llocal'))
user_parameters.append(optionify_params('broad', '--broad'))
user_parameters.append(optionify_params('broad_cutoff', '--broad-cutoff'))
user_parameters.append(optionify_params('cutoff_analysis', '--cutoff-analysis'))
user_parameters.append(optionify_params('call_summits', '--call-summits'))
user_parameters.append(optionify_params('fe_cutoff', '--fe-cutoff'))
user_parameters = ' '.join(user_parameters)
    
# Execute shell command.
shell(
    "("
    "macs2 callpeak "
    "{input_command} "
    "{output_command} "
    "{name_command} "
    "{user_parameters} "
    "{extra} "
    ") "
    "{log}"
)
