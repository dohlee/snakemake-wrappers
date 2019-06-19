rule bwameth_index:
    input:
        'hg38.fa'
    output:
        'hg38.fa.bwameth.c2t' 
    # NOTE: This rule does not require parameters.
    threads: 1  # NOTE: This rule will not use more than 1 core.
    log: 'logs/bwameth_index/_.log'
    benchmark: 'benchmarks/bwameth_index/_.log'
    wrapper:
        'http://dohlee-bio.info:9193/bwa-meth/index'

