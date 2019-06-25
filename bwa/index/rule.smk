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
        # Default: 'bwtsw',
        a = 'bwtsw',
        # Block size for the bwtsw algorithm (effective with -a bwtsw).
        # Default: False [10000000]
        b = False,
        # Index files named as <in.fasta>.64.* instead of <in.fasta>.*
        _6 = False,
    threads: 1
    log: 'logs/bwa_index/{genome}.log'
    benchmark: 'benchmarks/bwa_index/{genome}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bwa/index'
