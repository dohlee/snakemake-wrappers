from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'test.snpeff_annotated.vcf.gz'

rule fetch_vcf_data:
    input: HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/vcf/test.vcf.gz')
    output: 'test.vcf.gz'
    shell: 'mv {input} {output}'

rule snpeff:
    input:
        'test.vcf.gz'
    output:
        'test.snpeff_annotated.vcf.gz'
    params:
        # Required parameters.
        genome_version = 'hg38',
        # Optional parameters. Omit if unused.
        java_options = '-Xmx4g',
        # It true, there will be a significant speedup if there are a lot of samples.
        no_statistics = False,
        extra = ''
    threads: 1  # For now the wrapper only supports single-threaded execution.
    wrapper:
        'snpeff/ann'
