rule samtools_index:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.bam.bai'
    wrapper:
        'http://dohlee-bio.info:9193/samtools/index'
