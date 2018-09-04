from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'test.trimmed.fastq.gz'

rule trim_galore_se:
    input:
        HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/rrbs/se/test.fastq.gz')
    output:
        'test.trimmed.fastq.gz'
    threads: 1
    log:
        'logs/trim_galore/test.log'
    wrapper:
        'trim-galore'
