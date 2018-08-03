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

# Extract required arguments.
# Accession should start with 'SRR', 'ERR', or 'DRR'.
assert snakemake.output[0][:3] in ['SRR', 'ERR', 'DRR'], 'Accession should start with "SRR", "ERR", or "DRR".'
acc = snakemake.output[0].split('.')[0]

# Construct the path of sra file, according to the example sra path.
# ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR119/SRR1192353/SRR1192353.sra
ftp_path = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra' % (acc[:3], acc[:6], acc, acc)

# NOTE: Since it is hard to predict the output directory of `prefetch` command
# in sra-toolkit, I decided to just use `wget` to download sra file for now.
# Indeed we cannot enjoy the fast download speed of aspera client...
# TODO: Fix to use sra-toolkit `prefetch` command.
# Execute shell command.
shell(
    "("
    "wget "
    "-O {snakemake.output} "
    "{extra} "
    "{ftp_path} "
    ")"
    "{log}"
)
