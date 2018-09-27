rule sambamba_markdup:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.duplicates_marked.sorted.bam'
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/markdup'
