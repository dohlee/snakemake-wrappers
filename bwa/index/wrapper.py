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

# Assert input and output have been correctly given.
assert len(snakemake.input) == 1,
    'Please check your reference genome has been correctly given. It should be given as a single file.'

assert len(snakemake.output) == 5,
    'bwa-mem generates 5 outputs, *.amb, *.ann, *.bwt, *.pac, and *.sa. Please check your specified output.'

# Extract required inputs.
reference = snakemake.input[0]
prefix = path.join(path.dirname(reference), path.splitext(reference)[0])

algorithm = snakemake.params.algorithm.get('algorithm', 'bwtsw')
# Assert the algorithm is 'is' or 'bwtsw'.
assert algorithm in ['is', 'bwtsw'], 'Algorithm should be "is" or "bwtsw".'

# Execute shell command.
shell(
    "("
    "bwa index "
    "-p {prefix} "
    "-a {reference} "
    "{extra} "
    "{reference} "
    ") "
    "{log}"
)
