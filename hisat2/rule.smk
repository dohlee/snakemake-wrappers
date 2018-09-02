rule hisat2:
    input:
        # Required input.
        alignment = '{sample}.sorted.bam',
        annotation = 'reference/gencode.v28.gtf.gz'
    output:
        'result/{sample}/{sample}.htseq_count.result'
    threads: 1
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = '',
        # Mode to handle reads overlapping more thain one feature.
        # CANDIDATES: union, intersection-strict, intersection-nonempty
        # DEFAULT: union
        mode = 'union',
        # Whether the data is from a strand-specific assay.
        # CANDIDATES: yes, no, reverse
        # DEFAULT: yes
        # NOTE: Use reverse for the data generated recently.
        stranded = 'reverse'

    log: 'logs/htseq/count/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/hisat2'
