rule macs2_filterdup:
    input:
        '{sample}.sorted.bam'
    output:
        '{sample}.sorted.filterdup.bed'
    params:
        # MACS2 filterdup's behavior towards duplicate tags/pairs at the exact
        # same location.
        # 'auto': Calculate the maximum tags at the exact same location based on
        # binomial distribution.
        # integer value: Keep at most that much reads.
        # 'all': Keep all duplicates.
        keep_duplicate = 1,

        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        genome_size = 'hs',

        # Extra options.
        extra = ''
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat("benchmarks/macs2_filterdup/{sample}.tsv", 1)
    log: 'logs/macs2_filterdup/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/filterdup'

