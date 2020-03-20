rule sambamba_merge:
    input:
        ['run1.sorted.bam', 'run2.sorted.bam']
    output:
        '{sample}.sorted.bam'
    threads: 1
    log: 'logs/sambamba_merge/{sample}.log'
    benchmark: 'benchmarks/sambamba_merge/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/merge'
