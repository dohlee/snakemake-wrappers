rule sambamba_markdup:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.duplicates_marked.sorted.bam'
    threads: 1
    params:
        extra = '',
        # Remove duplicates instead of just marking them.
        # Default: False
        remove_duplicates = False,
        # Specify compression level of the resulting file (from 0 to 9).
        # Default: False
        compression_level = False,
        # Show progressbar in STDERR.
        # Default: False
        show_progress = False,
        # Specify directory for temporary files.
        # Default: False
        tmpdir = False,
        #
        # Performance tweaking parameters.
        #
        # Size of hash table for finding read pairs (default is 262144 reads);
        # will be rounded down to the nearest power of two;
        # should be (average coverage) * (insert size) for good performance.
        # Default: 262144
        hash_table_size = 262144,
        # Size of the overflow list where reads, thrown from the hash table,
        # get a second chance to meet their pairs (default is 200000 reads),
        # increasing the size reduces the number of temporary files created.
        # Default: 200000
        overflow_list_size = 200000,
        # Total amount of memory (in *megabytes*) used for sorting purposes;
        # the default is 2048, increasing it will reduce the number of created
        # temporary files and the time spent in the main thread.
        # Default: 2048
        sort_buffer_size = 2048,
        # Two buffers of BUFFER_SIZE *megabytes* each are used for reading and
        # writing BAM during the second pass.
        # Default: 128
        io_buffer_size = 128,
    log: 'logs/sambamba_markdup/{sample}.log'
    benchmark: 'benchmarks/sambamba_markdup/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/markdup'
