from os.path import join

rule bowtie_build:
    input:
        # Required input.
        # Reference genome fasta.
        # e.g.
        # reference = config['reference']['genome']['fasta']
        reference = 'REFERENCE_FASTA'
    output:
        # e.g.
        # reference_indices = config['reference']['genome']['bowtie_index_dir']
        index_dir = directory('DIRECTORY_FOR_YOUR_BOWTIE_INDEX')
    params:
        # Additional parameters go here.
        extra = '',
        # Build a colorspace index.
        color = False,
        # Disable automatic -p/--bmax/--dcv memory-fitting.
        noauto = False,
        # Use packed strings internally; slower, uses less mem.
        packed = False, 
        # Max bucket size for blockwise suffix-array builder.
        bmax = False,
        # Max bucket size as divisor of ref len.
        # Default: 4
        bmaxdivn = 4,
        # diff-cover period for blockwise
        # Default: 1024
        dcv = 1024,
        # disable diff-cover (algorithm becomes quadratic)
        nodc = False,
        # Don't build .3/.4.ebwt (packed reference) portion.
        noref = False,
        # Just build .3/.4.ebwt (packed reference) portion.
        justref = False,
        # SA is sampled every 2^offRate BWT chars.
        # Default: 5
        offrate = 5,
        # Number of chars consumed in initial lookup.
        # Default: 10
        ftabchars = 10,
        # convert Ns in reference to As.
        ntoa = False,
        # Seed for random number generator.
        seed = 0,
        # Verbose output (for debugging).
        quiet = False,
    threads: 4
    benchmark:
        repeat('benchmarks/bowtie_build/{genome}.tsv', 1)
    log: 'logs/bowtie_build/{genome}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bowtie/build'
