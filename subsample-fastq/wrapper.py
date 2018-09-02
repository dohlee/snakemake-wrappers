__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

def optionify_input(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.input[parameter])
        if param:
            return option + ' ' + str(snakemake.input[parameter])
        else:
            return ''
    except AttributeError:
        return ''

def optionify_params(parameter, option):
    """Return optionified parameter."""
    try:
        param = str(snakemake.params[parameter])
        if param:
            return option + ' ' + str(snakemake.params[parameter])
        else:
            return ''
    except AttributeError:
        return ''

# Extract required inputs.
reads = snakemake.input.reads
if isinstance(reads, str):
    reads = [reads]
assert len(reads) in [1, 2], "Input should be single-read or paired-end."

# Extract required outputs.
output = snakemake.output
gzip_required = all(o.endswith('.gz') for o in output)
ungzipped_output = [o[:-3] if o.endswith('.gz') else o for o in output]

# Extract parameters.
# Extract optional parameters.
k = snakemake.params.get('k', '')
assert k != '', 'Please provide the number of sampled reads (k).'

# Organize commands.
if all(r.endswith('.gz') for r in reads):
    if len(reads) == 1:
        cat = 'zcat %s' % reads[0]
    else:
        cat = 'paste <(zcat %s) <(zcat %s)' % (reads[0], reads[1])
else:
    if len(reads) == 1:
        cat = 'cat %s' % reads[0]
    else:
        cat = 'paste %s %s' % (reads[0], reads[1])

squish = 'awk \'{ printf("%s",$0); n++; if(n%4==0) { printf(\"\\n\");} else { printf(\"\\t\");} }\''
reservoir_sampling = 'awk -v k=%d \'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=NR<=k?NR:int(rand()*NR);if(s<k)R[s]=$0}END{for(i in R)print R[i]}\'' % k
if len(reads) == 1:
    print_out = 'awk -F\"\\t\" \'{print $1\"\\n\"$2\"\\n\"$3\"\\n\"$4 > \"%s\"}\'' % ungzipped_output[0]
    gzip = '&& gzip %s' % ungzipped_output[0] if gzip_required else ''
else:
    print_out = 'awk -F\"\\t\" \'{print $1\"\\n\"$3\"\\n\"$5\"\\n\"$7 > \"%s\"; print $2\"\\n\"$4\"\\n\"$6\"\\n\"$8 > \"%s\"}\'' % (ungzipped_output[0], ungzipped_output[1])
    gzip = '&& gzip %s && gzip %s' % (ungzipped_output[0], ungzipped_output[1]) if gzip_required else ''
# Execute shell command.
shell(
    "("
    "{cat} | "
    "{squish} | "
    "{reservoir_sampling} | "
    "{print_out} "
    "{gzip} "
    ") "
    "{log}"
)
