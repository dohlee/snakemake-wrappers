rule macs2_filterdup:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.filterdup.bed'
    params:
        # Extra options.
        extra = ''
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = 'hs',
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = False,
        # Pvalue cutoff for binomial distribution test.
        # Default: 1e-5
        pvalue = 1e-5,
        # MACS2 filterdup's behavior towards duplicate tags/pairs at the exact
        # same location.
        # 'auto': Calculate the maximum tags at the exact same location based on
        # binomial distribution.
        # integer value: Keep at most that much reads.
        # 'all': Keep all duplicates.
        keep_dup = 1,
        # Set verbose level.
        # 0: only show critical message.
        # 1: show additional warning message.
        # 2: show process information.
        # 3: show debug messages.
        # If you want to know where are the duplicate reads, use 3.
        # Default: 2
        verbose = 2,
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_filterdup/{sample}.tsv', 1)
    log: 'logs/macs2_filterdup/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/filterdup'

