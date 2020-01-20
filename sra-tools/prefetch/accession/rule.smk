rule prefetch_accession:
    output:
        temp('{sample}.sra')
    resources:
        network = 1
    wrapper:
        'http://dohlee-bio.info:9193/sra-tools/prefetch/accession'
