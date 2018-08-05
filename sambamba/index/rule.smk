rule samtools_sort:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.bam.bai'
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/index'
