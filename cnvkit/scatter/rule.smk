rule cnvkit_scatter:
    input:
        copy_ratio = ['{sample}.cnr'],
        segment = ['{sample}.cns'],
    output:
        '{sample}.cnv.pdf'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Plot segment lines in this color. value can be any string
        # accepted by matplotlib, e.g. 'red' or '#CC0000'
        segment_color = ''
        # Plot title. [Default: sample ID, from filename or -i]
        title = '',
    threads: 1
    log: 'logs/cnvkit/scatter/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/scatter'
