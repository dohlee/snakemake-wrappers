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

def is_defined_by_user(*params):
    extra = snakemake.params.get('extra', '')
    for param in params:
        if param in extra:
            return True
    return False

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
user_parameters.append(optionify_params('seed', '--seed'))
user_parameters.append(optionify_params('skip', '--skip'))
user_parameters.append(optionify_params('qupto', '--qupto'))
user_parameters.append(optionify_params('trim5', '--trim5'))
user_parameters.append(optionify_params('trim3', '--trim3'))
user_parameters.append(optionify_params('phred33_quals', '--phred33-quals'))
user_parameters.append(optionify_params('phred64_quals', '--phred64-quals'))
user_parameters.append(optionify_params('solexa_quals', '--solexa-quals'))
user_parameters.append(optionify_params('solexa13_quals', '--solexa1.3-quals'))
user_parameters.append(optionify_params('integer_quals', '--integer-quals'))
user_parameters.append(optionify_params('large_index', '--large-index'))
user_parameters.append(optionify_params('v', '-v'))
user_parameters.append(optionify_params('seedmms', '--seedmms'))
user_parameters.append(optionify_params('maqerr', '--maqerr'))
user_parameters.append(optionify_params('seedlen', '--seedlen'))
user_parameters.append(optionify_params('nomaqround', '--nomaqround'))
user_parameters.append(optionify_params('minins', '--minins'))
user_parameters.append(optionify_params('maxins', '--maxins'))
user_parameters.append(optionify_params('fr', '--fr'))
user_parameters.append(optionify_params('rf', '--rf'))
user_parameters.append(optionify_params('ff', '--ff'))
user_parameters.append(optionify_params('nofw', '--nofw'))
user_parameters.append(optionify_params('norc', '--norc'))
user_parameters.append(optionify_params('maxbts', '--maxbts'))
user_parameters.append(optionify_params('pairtries', '--pairtries'))
user_parameters.append(optionify_params('tryhard', '--tryhard'))
user_parameters.append(optionify_params('chunkmbs', '--chunkmbs'))
user_parameters.append(optionify_params('reads_per_batch', '--reads-per-batch'))
user_parameters.append(optionify_params('k', '-k'))
user_parameters.append(optionify_params('all', '--all'))
user_parameters.append(optionify_params('m', '-m'))
user_parameters.append(optionify_params('M', '-M'))
user_parameters.append(optionify_params('best', '--best'))
user_parameters.append(optionify_params('strata', '--strata'))
user_parameters.append(optionify_params('time', '--time'))
user_parameters.append(optionify_params('offbase', '--offbase'))
user_parameters.append(optionify_params('quiet', '--quiet'))
user_parameters.append(optionify_params('refidx', '--refidx'))
user_parameters.append(optionify_params('al', '--al'))
user_parameters.append(optionify_params('un', '--un'))
user_parameters.append(optionify_params('no_unal', '--no-unal'))
user_parameters.append(optionify_params('max', '--max'))
user_parameters.append(optionify_params('suppress', '--suppress'))
user_parameters.append(optionify_params('fullref', '--fullref'))
user_parameters.append(optionify_params('snpphred', '--snpphred'))
user_parameters.append(optionify_params('snpfrac', '--snpfrac'))
user_parameters.append(optionify_params('col_cseq', '--col-cseq'))
user_parameters.append(optionify_params('col_cqual', '--col-cqual'))
user_parameters.append(optionify_params('col_keepends', '--cool-keepends'))
user_parameters.append(optionify_params('sam', '--sam'))
user_parameters.append(optionify_params('mapq', '--mapq'))
user_parameters.append(optionify_params('sam_nohead', '--sam-nohead'))
user_parameters.append(optionify_params('sam_nosq', '--sam-nosq'))
user_parameters.append(optionify_params('sam_rg', '--sam-RG'))
user_parameters.append(optionify_params('offrate', '--offrate'))
user_parameters.append(optionify_params('mm', '--mm'))
user_parameters.append(optionify_params('shmem', '--shmem'))
user_parameters = ' '.join([p for p in user_parameters if p != ' '])

# Extract required inputs.
reads = snakemake.input.reads
if len(reads) == 2:
    read_command = '-1 %s -2 %s' % (reads[0], reads[1])
else:
    read_command = '%s' % (reads[0])

index_dir = snakemake.input.index_dir
prefix = get_common_prefixes([f for f in listdir(index_dir) if not f.startswith('.')])
index_command = path.join(index_dir, prefix)

# Extract required outputs.
output = snakemake.output[0]
if output.endswith('.sorted.bam'):
    if 'mapq_cutoff' in snakemake.params:
        postprocess_command = '| samtools view -bS -q %d - | samtools sort > %s' % (snakemake.params.mapq_cutoff, output)
    else:
        postprocess_command = '| samtools view -bS - | samtools sort > %s' % (output)
else:
    if 'mapq_cutoff' in snakemake.params:
        postprocess_command = '| samtools view -bS -q %d - > %s' % (snakemake.params.mapq_cutoff, output)
    else:
        postprocess_command = '| samtools view -bS - > %s' % (output)

# Execute shell command.
shell(
    "("
    "bowtie "
    "{index_command} "
    "{read_command} "
    "{extra} "
    "{user_parameters} "
    "--threads {snakemake.threads} "
    "{postprocess_command}) "
    "{log}"
)
