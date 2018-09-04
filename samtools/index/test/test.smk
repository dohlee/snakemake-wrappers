from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'test.sorted.bam.bai'

rule samtools_sort:
    input:
        HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/bam/test.bam')
    output:
        'test.sorted.bam'
    wrapper:
        'samtools/sort'

rule samtools_index:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.bam.bai'
    wrapper:
        'samtools/index'
