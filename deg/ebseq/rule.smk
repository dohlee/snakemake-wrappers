rule ebseq:
    input:
        # Required.
        data = 'GENE_EXPRESSION_RAW_COUNT_TABLE',
        condition = 'CONDITION_TABLE',
    threads: 1
    output:
        # Required.
        deg_list = 'result/ebseq_deg.list',
        fold_change = 'result/ebseq_foldchange.tsv',
        result = 'ebseq_result.tsv',
    params:
        cutoff = 0.05,
        verbose = False,
    wrapper:
        'http://dohlee-bio.info:9193/deg/ebseq'
