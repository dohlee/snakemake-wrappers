__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2023, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract required inputs.
pdb_id = os.path.splitext(os.path.basename(snakemake.output))[0]

# Extract required outputs.
output = snakemake.output

# Execute shell command.
shell(
    "("
    "wget "
    "https://files.rcsb.org/download/{pdb_id}.pdb "
    "-O {output}"
    ") "
    "{log}"
)
