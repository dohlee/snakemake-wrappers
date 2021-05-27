__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
import os
# NOTE: Since it is hard to predict the output directory of `prefetch` command
# in sra-toolkit, I decided to just use `wget` to download sra file for now.
# Indeed we cannot enjoy the fast download speed of aspera client...
# TODO: Fix to use sra-toolkit `prefetch` command.

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
# Accession should start with 'SRR', 'ERR', or 'DRR'.
assert snakemake.output[0][:3] in ['SRR', 'ERR', 'DRR'], 'Accession should start with "SRR", "ERR", or "DRR".'
acc = os.path.basename(snakemake.output[0]).split('.')[0]
out = snakemake.output

command_template = 'wget https://sra-download{region}.ncbi.nlm.nih.gov/sos{sos}/sra-pub-run-{run}/' + str(acc) + '/' + str(acc) + '.{i} -qO ' + str(out)

wget_commands = []
for region in ['b.be-md', '.st-va']:
    for sos in [1, 2]:
        for run in range(1, 11):
            for i in [1, 2, 3]:
                wget_command = command_template.format(region=region, sos=sos, run=run, i=i)
                wget_commands.append(wget_command)

command = ' || '.join(wget_commands)

# Execute shell command.
shell(
    "("
    "{command}"
    ")"
    "{log}"
)
