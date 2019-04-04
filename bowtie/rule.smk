rule bowtie:
    input:
        index_dir = directory('DIRECTORY_FOR_YOUR_BOWTIE_INDEX'),
        # For single reads,
        #reads = [
        #    '{sample}.fastq.gz',
        #],
        # For paired-end reads,
        reads = [
            '{sample}.read1.fastq.gz',
            '{sample}.read2.fastq.gz',
        ]
    output:
        # It automatically sorts the output bam file if its file name ends with '.sorted.bam',
        # e.g.
        # '{sample}.sorted.bam'
        bam = '{sample}.bam'
    params:
        # Additional parameters go here.
        extra = '',
        # Discard mapped reads having mapping quality (MAPQ) below this value.
        mapq_cutoff = 10,
    threads: 4
    benchmark:
        repeat("benchmarks/bowtie/{sample}.tsv", 1)
    log: 'logs/bowtie/{sample}.log'
    wrapper:
    'http://dohlee-bio.info:9193/bowtie'

