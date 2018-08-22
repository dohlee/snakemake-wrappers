rule ebseq:
    input:
        # Required.
        data = '',
        condition = '',
    threads: 1
    output:
        # Required.
        deg_list = '',
        result = '',
    params:
        cutoff = 0.05,
        verbose = False,
    wrapper:
        'http://dohlee-bio.info:9193/deg/deseq2'
