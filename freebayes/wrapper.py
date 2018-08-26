__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required inputs.
bams = snakemake.input.bam
bams = [snakemake.input] if isinstance(bams, str) else bams
bams_option = ' '.join(['--bam ' + f for f in bams])

reference = snakemake.input.reference

# Extract required outputs.
output = snakemake.output[0]

# If user wants bcf file as an output, pipe the output through bcftools.
pipe_command = ''
if output.endswith('.bcf'):
    pipe_command += '| bcftools view -Ob -'

chunksize = snakemake.params.get('chunksize', 100000)

if snakemake.threads == 1:
    freebayes = 'freebayes'
else:
    freebayes = 'freebayes-parallel <(fasta_generate_regions.py %s.fai %d %d)' % (reference, chunksize, snakemake.threads)

# Execute shell command.
shell(
    "("
    "{freebayes} "
    "{extra} "
    "-f {reference} "
    "{bams_option} "
    ") "
    "{pipe_command} > "
    "{output} "
    "{log}"
)
