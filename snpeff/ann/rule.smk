rule snpeff:
    input:
        '{sample}.vcf.gz'
    output:
        '{sample}.snpeff_annotated.vcf.gz'
    params:
        # Required parameters.
        genome_version = 'hg38',
        # Optional parameters. Omit if unused.
        java_options = '-Xmx4g',
        # It true, there will be a significant speedup if there are a lot of samples.
        no_statistics = False,
        extra = ''
    threads: 1  # For now the wrapper only supports single-threaded execution.
    log: 'logs/snpeff/ann/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/snpeff/ann'
