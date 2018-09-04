__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
reference_dir = snakemake.input

# Execute shell command.
shell(
    "("
    "bismark_genome_preparation "
    "{extra} "
    "{reference_dir} "
    ")"
    "{log}"
)
