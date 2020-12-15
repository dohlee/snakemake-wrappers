rule to_tdf:
    input:
        '{sample}.bdg'
    output:
        '{sample}.bdg.tdf'
    params:
        genome_version = 'hg38'
    threads: 1
    log: 'logs/igvtools/to_tdf/{sample}.log'
    benchmark: repeat('benchmarks/igvtools/to_tdf/{sample}.benchmark', 1)
    wrapper:
        'http://dohlee-bio.info:9193/igvtools/totdf'
