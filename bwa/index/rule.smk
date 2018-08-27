rule bwa_index:
    input:
        # Required input. Reference genome fasta file.
        '{genome}.fasta'
    output:
        # Required output. BWA-indexed reference genome files.
        '{genome}.amb',
        '{genome}.ann',
        '{genome}.bwt',
        '{genome}.pac',
        '{genome}.sa'
    params:
        extra = '',
        # Note that the default algorithm for this wrapper is 'bwtsw'.
        # The other option is 'is', but please be warned that this algorithm doesn't work
        # with genomes longer than 2GB.
        algorithm = 'bwtsw'
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/bwa/index'
