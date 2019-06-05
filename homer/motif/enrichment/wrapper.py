__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from os import dirname
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
user_parameters.append(optionify_params('mask', '-mask'))
user_parameters.append(optionify_params('bg', '-bg'))
user_parameters.append(optionify_params('len_', '-len'))
user_parameters.append(optionify_params('size', '-size'))
user_parameters.append(optionify_params('S', '-S'))
user_parameters.append(optionify_params('mis', '-mis'))
user_parameters.append(optionify_params('norevopp', '-norevopp'))
user_parameters.append(optionify_params('nomotif', '-nomotif'))
user_parameters.append(optionify_params('rna', '-rna'))
user_parameters.append(optionify_params('mset', '-mset'))
user_parameters.append(optionify_params('basic', '-basic'))
user_parameters.append(optionify_params('bits', '-bits'))
user_parameters.append(optionify_params('nocheck', '-nocheck'))
user_parameters.append(optionify_params('mcheck', '-mcheck'))
user_parameters.append(optionify_params('float_', '-float'))
user_parameters.append(optionify_params('noknown', '-noknown'))
user_parameters.append(optionify_params('mknown', '-mknown'))
user_parameters.append(optionify_params('nofacts', '-nofacts'))
user_parameters.append(optionify_params('seqlogo', '-seqlogo'))
user_parameters.append(optionify_params('gc', '-gc'))
user_parameters.append(optionify_params('cpg', '-cpg'))
user_parameters.append(optionify_params('noweight', '-noweight'))
user_parameters.append(optionify_params('h', '-h'))
user_parameters.append(optionify_params('N', '-N'))
user_parameters.append(optionify_params('local', '-local'))
user_parameters.append(optionify_params('redundant', '-redundant'))
user_parameters.append(optionify_params('maxN', '-maxN'))
user_parameters.append(optionify_params('maskMotif', '-maskMotif'))
user_parameters.append(optionify_params('opt', '-opt'))
user_parameters.append(optionify_params('rand', '-rand'))
user_parameters.append(optionify_params('ref', '-ref'))
user_parameters.append(optionify_params('oligo', '-oligo'))
user_parameters.append(optionify_params('dumpFasta', '-dumpFasta'))
user_parameters.append(optionify_params('preparse', '-preparse'))
user_parameters.append(optionify_params('preparsedDir', '-preparsedDir'))
user_parameters.append(optionify_params('keepFiles', '-keepFiles'))
user_parameters.append(optionify_params('fdr', '-fdr'))
user_parameters.append(optionify_params('homer2', '-homer2'))
user_parameters.append(optionify_params('nlen', '-nlen'))
user_parameters.append(optionify_params('nmax', '-nmax'))
user_parameters.append(optionify_params('neutral', '-neutral'))
user_parameters.append(optionify_params('olen', '-olen'))
user_parameters.append(optionify_params('e', '-e'))
user_parameters.append(optionify_params('cache', '-cache'))
user_parameters.append(optionify_params('quickMask', '-quickMask'))
user_parameters.append(optionify_params('minlp', '-minlp'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required inputs.
bed = snakemake.input.bed
genome = snakemake.input.genome

# Extract required outputs.
homer_results = snakemake.output.homer_results
known_results = snakemake.output.known_results

outdir = dirname(homer_results)

# Execute shell command.
shell(
    "("
    "findMotifsGenome.pl "
    "{bed} "
    "{genome} "
    "{outdir} "
    "{user_parameters} "
    "-p {snakemake.threads} "
    ")"
    "{log}"
)
