rule parallel_fastq_dump_single:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        '{sample}.fastq.gz'
    params:
        extra = '--tmpdir .'
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'

rule parallel_fastq_dump_paired:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        '{sample}.read1.fastq.gz',
        '{sample}.read2.fastq.gz',
        temp('{sample}_pass.fastq.gz')
    params:
        # Optional parameters. Omit if unused.
        extra = '--tmpdir .'
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'
