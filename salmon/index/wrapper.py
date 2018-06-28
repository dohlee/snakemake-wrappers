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
fasta = snakemake.input[0]
index = snakemake.output[0]
threads = snakemake.threads

# Execute shell command.
shell(
    "("
    "salmon index "
    "--transcripts {fasta} "
    "--index {index} "
    "--threads {threads} "
    "{extra} "
    ")"
    "{log}"
)
