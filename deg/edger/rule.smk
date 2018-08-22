rule edgeR:
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
        dispersion = 'common',  # common, trended, tagwise.
        verbose = False,
    wrapper:
        'http://dohlee-bio.info:9193/deg/edger'
