__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
output = snakemake.output[0]
BASE_URL = 'https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/annotation/'

if 'hg38' in output:
    url = BASE_URL + 'Homo_sapiens.hg38.annotation.gtf'
elif 'GRCh38' in output:
    url = BASE_URL + 'Homo_sapiens.GRCh38.annotation.gtf'
else:
    raise ValueError('Output annotation file name should include "hg38" or "GRCh38".')

# Execute shell command.
shell(
    "("
    "wget {url} -qO {output}"
    ")"
)
