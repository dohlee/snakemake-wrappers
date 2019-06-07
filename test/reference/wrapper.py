__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
output = snakemake.output[0]
if output.endswith('.gz'):
    wget_output = output
    gunzip_command = ':'
else:
    wget_output = output + '.gz'
    gunzip_command = 'gunzip %s' % wget_output

BASE_URL = 'https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/reference/'
if 'GRCh38' in output:
    url = BASE_URL + 'Homo_sapiens.GRCh38.genome.chr20.fasta.gz'
elif 'hg38' in output:
    url = BASE_URL + 'Homo_sapiens.hg38.genome.chr20.fasta.gz'
elif 'mm10' in output:
    url = BASE_URL + 'Mus_musculus.mm10.genome.chr19.fasta.gz'
else:
    raise ValueError('Please specify "GRCh38" or "hg38" or "mm10" in your output file.')

# Execute shell command.
shell(
    "("
    "wget {url} -qO {wget_output} && "
    "{gunzip_command}"
    ")"
)
