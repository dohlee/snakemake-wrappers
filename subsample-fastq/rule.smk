rule subsample_fastq:
    input:
        # Required input.
        reads = ['data/{sample}/{sample}.fastq.read1.gz', 'data/{sample}/{sample}.fastq.read2.gz']
    output:
        'data/{sample}/{sample}.subsampled.read1.fastq.gz',
        'data/{sample}/{sample}.subsampled.read2.fastq.gz'
    threads: 1  # No more than 1 threads will be used.
    params:
        # Required parameters.
        k = 10000  # Number of sampled reads.
    log: 'logs/fastq_subsampling/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/subsample-fastq'
