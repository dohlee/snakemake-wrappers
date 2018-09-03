rule sambamba_markdup:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.duplicates_marked.bam'
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/markdup'
