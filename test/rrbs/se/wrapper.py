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

url = 'https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rrbs/se/test.fastq.gz'

# Execute shell command.
shell(
    "("
    "wget {url} -qO {wget_output} && "
    "{gunzip_command}"
    ")"
)
