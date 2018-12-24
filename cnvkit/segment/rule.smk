rule cnvkit_segment:
    input:
        copy_ratios = '{tumor}.cnr',
    output:
        '{tumor}.cns',
    params:
        extra = '',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/segment'
