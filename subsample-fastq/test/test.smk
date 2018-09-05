from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'test.subsampled.fastq.gz', 'test1.subsampled.fastq'

rule download_data:
    input: HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/rna-seq/se/test.fastq.gz')
    output: 'test.fastq.gz'
    shell: 'mv {input} {output}'

rule gunzip:
    input: 'test.fastq.gz'
    output: 'test1.fastq'
    shell: 'gunzip -c {input} > {output}'

rule subsample_fastq_gzipped_se:
    input:
        # Required input.
        reads = ['test.fastq.gz']
    output:
        'test.subsampled.fastq.gz',
    threads: 1  # No more than 1 threads will be used.
    params:
        # Required parameters.
        k = 1000  # Number of sampled reads.
    log: 'logs/fastq_subsampling/test.log'
    wrapper:
        'subsample-fastq'

rule subsample_fastq_se:
    input:
        # Required input.
        reads = ['test1.fastq']
    output:
        'test1.subsampled.fastq',
    threads: 1  # No more than 1 threads will be used.
    params:
        # Required parameters.
        k = 1000  # Number of sampled reads.
    log: 'logs/fastq_subsampling/test.log'
    wrapper:
        'subsample-fastq'
