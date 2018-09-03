rule deseq2:
    input:
        # Required.
        data = 'GENE_EXPRESSION_RAW_COUNT_TABLE',
        condition = 'CONDITION_TABLE',
    threads: 1
    output:
        # Required.
        deg_list = 'result/deseq2_deg.list',
        result = 'result/deseq2_result.tsv',
    params:
        cutoff = 0.05,
        verbose = False,
    wrapper:
        'http://dohlee-bio.info:9193/deg/deseq2'
