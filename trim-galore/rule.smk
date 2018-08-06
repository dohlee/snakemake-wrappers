rule trim_galore_se:
    input:
        '{sample}.fastq.gz'
    output:
        '{sample}.trimmed.fastq.gz'
    threads: 1
    log:
        'logs/trim_galore/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/trim-galore'

rule trim_galore_pe:
    input:
        '{sample}.read1.fastq.gz',
        '{sample}.read2.fastq.gz'
    output:
        '{sample}.read1.trimmed.fastq.gz',
        '{sample}.read2.trimmed.fastq.gz'
    threads: 1
    log:
        'logs/trim_galore/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/trim-galore'
