rule sambamba_index:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.bam.bai'
    threads: 1
    log: 'logs/sambamba_index/{sample}.log'
    benchmark: 'benchmarks/sambamba_index/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/index'
