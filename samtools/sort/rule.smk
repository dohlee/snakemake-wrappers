rule samtools_sort:
    input:
        '{sample}.bam'
    output:
        '{sample}.sorted.bam'
    params:
        extra = ''
        # Set compression level.
        # '-l INT '
        # Set maximum memory per thread; suffix K/M/G/ recognized [768M].
        # '-m INT '
        # Sort by read name.
        # '-n '
        # Sort by value of TAG.
        # '-t TAG '
    wrapper:
        'http://dohlee-bio.info:9193/samtools/sort'
