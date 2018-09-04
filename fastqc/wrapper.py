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
read_command = ' '.join(snakemake.input)

html = snakemake.output.html
zipped = snakemake.output.zip
output_directory = path.dirname(html)

# Execute shell command.
shell(
    "("
    "fastqc "
    "-o {output_directory} "
    "-t {snakemake.threads} "
    "{extra} "
    "{read_command} "
    ") "
    "{log}"
)
