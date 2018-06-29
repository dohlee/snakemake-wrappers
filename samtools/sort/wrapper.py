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
input = snakemake.input[0]
output = snakemake.output[0]
prefix = path.splitext(output)[0]

# Execute shell command.
shell(
    "("
    "samtools sort "
    "{extra} "
    "-@ {snakemake.threads} "
    "-o {output} "
    "-T {prefix} "
    "{input}"
    ")"
    "{log}"
)
