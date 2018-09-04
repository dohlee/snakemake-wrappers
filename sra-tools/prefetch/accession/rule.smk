rule prefetch_accession:
    output:
        'SRR******.sra'
    wrapper:
        'http://dohlee-bio.info:9193/sra-tools/prefetch/accession'
