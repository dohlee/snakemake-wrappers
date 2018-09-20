rule samtools_faidx:
    input:
        '{file}.fasta'
    output:
        '{file}.fasta.fai'
    wrapper:
        'http://dohlee-bio.info:9194/samtools/faidx'
