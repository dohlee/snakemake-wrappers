rule sambamba_sort:
    input:
        '{sample}.bam'
    output:
        '{sample}.sorted.bam'
    params:
        extra = '',
        # Sets an upper bound for used memory. However, this is very approximate.
        # Default memory limit is 512 MiB. Increasing it will allow to make chunk
        # sizes larger and also reduce amount of I/O seeks thus improving the
        # overall performance. LIMIT must be a number with an optional suffix
        # specifying unit of measurement.
        # The following endings are recongnized: K, KiB, KB, M, MiB, MB, G, GiB, GB
        # Default: False
        memory_limit = False,
        # Use TMPDIR to output sorted chunks. Default behavior is to use system 
        # temporary directory.
        # Default: False
        tmpdir = False,
        # Sort by read name instead of doing coordinate sort.
        # Default: False
        sort_by_name = False,
        # Compression level to use for sorted BAM, from 0 (known as uncompressed
        # BAM in samtools) to 9.
        # Default: False
        compression_level = False,
        # Write sorted chuks as uncompressed BAM. Default dehaviour is to write them
        # with compression level 1, because that reduces time spent on I/O, but in
        # some cases using this option can give you a better speed. Note, however,
        # that the disk space needed for sorgin will typically be 3-4 times more
        # than without enabling this option.
        # Default: False
        uncompressed_chunks = False,
        # Show wget-like progressbar in STDERR (in fact, two of them one after another,
        # first one for sorting, and then another one for merging.)
        # Default: False
        show_progress = False,
    threads: 1
    log: 'logs/sambamba_sort/{sample}.log'
    benchmark: 'benchmarks/sambamba_sort/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/sort'
