__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Define exception classes.
class RuleInputException(Exception):
    pass

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

# Extract parameters.
extra = snakemake.params.get('extra', '')
user_parameters = []
user_parameters.append(optionify_params('gff3_rna_patterns', '--gff3-RNA-patterns'))
user_parameters.append(optionify_params('trusted_sources', '--trusted-sources'))
user_parameters.append(optionify_params('transcript_to_gene_map', '--transcipt-to-gene-map'))
user_parameters.append(optionify_params('allele_to_gene_map', '--allele-to-gene-map'))
user_parameters.append(optionify_params('quiet', '--quiet'))
user_parameters = ' '.join([p for p in user_parameters if p != ''])

# Extract required arguments.
fasta = snakemake.input.fasta
gtf = snakemake.input.get('gtf', None)
gff3 = snakemake.input.get('gff3', None)

if (gtf is not None) and (gff3 is not None):
    raise RuleInputException('You cannot provide both GTF and GFF3 files.')
annotation_option = ''
if gtf is not None:
    annotation_option = '--gtf %s' % gtf
elif gff3 is not None:
    annotation_option = '--gff3 %s' % gff3

output_prefix = snakemake.output[0].rstrip('.transcripts.fa')
threads = snakemake.threads

# Execute shell command.
shell(
    "("
    "rsem-prepare-reference "
    "{annotation_option} "
    "{extra} "
    "{user_parameters} "
    "{fasta} "
    "{output_prefix} "
    ")"
    "{log}"
)
