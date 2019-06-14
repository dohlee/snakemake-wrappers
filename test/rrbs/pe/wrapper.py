__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
read1 = snakemake.output[0]
read2 = snakemake.output[1]

if read1.endswith('.gz'):
    read1_wget_output = read1
    read1_gunzip_command = ':'
else:
    read1_wget_output = read1 + '.gz'
    read1_gunzip_command = 'gunzip %s' % wget_output

if read2.endswith('.gz'):
    read2_wget_output = read2
    read2_gunzip_command = ':'
else:
    read2_wget_output = read2 + '.gz'
    read2_gunzip_command = 'gunzip %s' % wget_output

# Execute shell command.
shell(
    "("
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rrbs/pe/test.read1.fastq.gz -qO {read1_wget_output} && "
    "{read1_gunzip_command} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rrbs/pe/test.read2.fastq.gz -qO {read2_wget_output} && "
    "{read2_gunzip_command}"
    ")"
)
