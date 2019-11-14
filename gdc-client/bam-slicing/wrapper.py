__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract required arguments.
token_file = snakemake.input.token
with open(token_file) as inFile:
    token = inFile.readline().strip()

out = snakemake.output[0]

parameters = []
if 'region' in snakemake.params:
    for r in snakemake.params['region'].split(','):
        parameters.append('region=' + r)
if 'gencode' in snakemake.params:
    for g in snakemake.params['gencode'].split(','):
        parameters.append('gencode=' + g)

if len(parameters) != 0:
    raise ValueError('Please specify at least one parameters for BAM slicing.')
parameters = '&'.join(parameters)

# Execute shell command.
shell(
    '('
    'curl --header "X-Auth-Token: {token}" '
    '\'https://api.gdc.cancer.gov/slicing/view/{snakemake.params.uuid}?{parameters}\' '
    '--output {out}'
    ') '
    '{log}'
)
