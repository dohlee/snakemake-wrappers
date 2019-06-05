rule reference_transcriptome:
    output:
        'GRCh38_transcriptome_chr20.fa.gz'
    wrapper:
        'http://dohlee-bio.info/test/reference/transcriptome'
