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
user_parameters.append(('gsize', '--gsize'))
user_parameters.append(('keep_dup', '--keep-dup'))
user_parameters.append(('buffer_size', '--buffer-size'))
user_parameters.append(('name', '--name'))
user_parameters.append(('bdg', '--bdg'))
user_parameters.append(('trackline', '--trackline'))
user_parameters.append(('SPMR', '--SPMR'))
user_parameters.append(('tsize', '--tsize'))
user_parameters.append(('bw', '--bw'))
user_parameters.append(('mfold', '--mfold'))
user_parameters.append(('fix_bimodal', '--fix-bimodal'))
user_parameters.append(('nomodel', '--nomodel'))
user_parameters.append(('shift', '--shift'))
user_parameters.append(('extsize', '--extsize'))
user_parameters.append(('qvalue', '--qvalue'))
user_parameters.append(('pvalue', '--pvalue'))
user_parameters.append(('to_large', '--to-large'))
user_parameters.append(('ratio', '--ratio'))
user_parameters.append(('down_sample', '--down-sample'))
user_parameters.append(('seed', '--seed'))
user_parameters.append(('tempdir', '--tempdir'))
user_parameters.append(('nolambda', '--nolambda'))
user_parameters.append(('slocal', '--slocal'))
user_parameters.append(('llocal', '--llocal'))
user_parameters.append(('broad', '--broad'))
user_parameters.append(('broad_cutoff', '--broad-cutoff'))
user_parameters.append(('cutoff_analysis', '--cutoff-analysis'))
user_parameters.append(('call_summits', '--call-summits'))
user_parameters.append(('fe_cutoff', '--fe-cutoff'))
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
