rule edgeR:
    input:
        # Required.
        data = 'GENE_EXPRESSION_RAW_COUNT_TABLE',
        condition = 'CONDITION_TABLE',
    threads: 1
    output:
        # Required.
        deg_list = 'result/edger_deg.list',
        result = 'result/edger_result.tsv',
    params:
        cutoff = 0.05,
        dispersion = 'common',  # common, trended, tagwise.
        verbose = False,
    wrapper:
        'http://dohlee-bio.info:9193/deg/edger'
