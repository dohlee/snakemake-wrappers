rule snpsift_annotate:
    input:
        vcf = '{sample}.vcf.gz',
        vcf_index = '{sample}.vcf.gz.tbi',
        db = 'dbsnp151.vcf.gz',
    output:
        '{sample}.snpsift.vcf.gz'
    params:
        # Extra parameters.
        extra = '',
        java_options = '',
        # Required parameters.
        # Only annotate ID field (do not add INFO field).
        # Default: True
        id = True,
        # Annotate using a list of info fields (list is a comma separated list of fields).
        # Default: ALL
        info = False,
        # VCF database is tabix-indexed.
        # Default: False
        tabix = False,
    threads: 1  # For now the wrapper only supports single-threaded execution.
    log: 'logs/snpsift/annotate/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/snpsift/annotate'
