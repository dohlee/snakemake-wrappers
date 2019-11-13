rule sambamba_flagstat:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.bam.flagstat'
    threads: 1
    benchmark: 'benchmarks/sambamba_flagstat/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/flagstat'
