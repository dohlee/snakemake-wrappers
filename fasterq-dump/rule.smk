rule fasterq_dump_single:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        '{sample}.fastq.gz'
    params:
        extra = ''
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/fasterq-dump'

rule fasterq_dump_paired:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        '{sample}.read1.fastq.gz',
        '{sample}.read2.fastq.gz',
    params:
        # Optional parameters. Omit if unused.
        extra = ''
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/fasterq-dump'
