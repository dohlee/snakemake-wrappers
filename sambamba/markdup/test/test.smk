from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'test.sorted.duplicates_marked.bam'

rule sambamba_sort:
    input:
        HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/bam/test.bam')
    output:
        'test.sorted.bam'
    wrapper:
        'sambamba/sort'

rule sambamba_markdup:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.duplicates_marked.bam'
    threads: 1
    wrapper:
        'sambamba/markdup'
