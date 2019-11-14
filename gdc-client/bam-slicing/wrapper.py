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
is_region_specified = snakemake.params.get('region', None) is not None
is_gencode_specified = snakemake.params.get('gencode', None) is not None

print(snakemake.params.get('region', None))
print(snakemake.params['region'])
if is_region_specified:
    for r in snakemake.params['region'].split(','):
        parameters.append('region=' + r)
if is_gencode_specified:
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
