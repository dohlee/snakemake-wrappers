rule htseq_count:
    input:
        # Required input.
        alignment = '{sample}.sorted.bam',
        annotation = 'reference/gencode.v28.gtf.gz',
    output:
        'result/{sample}/{sample}.htseq_count.result',
    threads: 1
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = '',
        # Sorting order of alignment file.
        # Option: 'pos' or 'name'
        # Default: name
        order = 'name',
        # When alignment file is paired end sorted by position, allow only
        # so many reads to stay in memory until the mates are found.
        # (raising this number will user more memory). Has no effect for
        # single end or paired end sorted by name.
        max_reads_in_buffer = False,
        # Whether the data is from a strand-specific assay.
        # 'reverse' means 'yes' with reversed strand interpretation.
        # Option: 'yes', 'no', or 'reverse'
        # Default: yes
        stranded = 'reverse',
        # Skip all reads with alignment quality lower than the given
        # minimum value
        # Default: 10
        minaqual = 10,
        # Feature type (3rd column in GFF file) to be used, all
        # features of other type are ignored.
        # Default, suitable for Ensembl GTF files: exon
        type = False,
        # GFF attribute to be used as feature ID.
        # Default, suitable for Ensembl GTF files: gene_id
        idattr = False,
        # Additional feature attributes.
        # Use multiple times for each different attribute.
        # Default, none, suitable for Ensembl GTF files: gene_name
        additional_attr = False,
        # Mode to handle reads overlapping more than one feature.
        # Option: union, intersection-strict, intersection-nonempty
        # Default: union
        mode = 'union',
        # Whether to score reads that are not uniquely aligned or ambiguously assigned
        # to features.
        # Option: none, all
        nonunique = False,
        # Whether to score secondary alignments (0x100 flag)
        # Option: score, ignore
        secondary_alignments = False,
        # Whether to score supplementary alignments (0x800 flag)
        supplementary_alignments = False,
        # Write out all SAM alignment records into SAM files (one per input file needed),
        # annotating each line with its feature assignment (as an optional field with tag 'XF')
        samout = False,
        # Suppress progress report.
        quiet = False,
    log: 'logs/htseq/count/{sample}.log'
    benchmark: 'benchmarks/htseq/count/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/htseq/count'
