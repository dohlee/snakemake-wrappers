rule samtools_sort:
    input:
        '{sample}.bam'
    output:
        '{sample}.sorted.bam'
    threads: 1
    params:
        extra = ''
        # Approximate total memory limit for all threads [2GB].
        # '-m INT '
        # Directory for storing intermediate files.
        # '--tmpdir=TMPDIR '
        # Sort by read name.
        # '-n '
        # Compression level for sorted BAM, from 0 to 9.
        # '--compression-leve=COMPRESSION_LEVEL '
        # Keep only reads that satisfy FILTER.
        # '--filter=FILTER '
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/sort'
