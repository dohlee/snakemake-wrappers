rule fastqc:
    input:
        '{prefix}.fastq.gz'
    output:
        html='{prefix}_fastqc.html',
        zip='{prefix}_fastqc.zip',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/fastqc'
