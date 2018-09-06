rule feature_counts:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        alignment = ['{sample}.bam'],
        annotation = 'annotation/gencode.v28.gtf.gz'
    output:
        # Required output.
        '{sample}.feature_counts.result'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Format of provided annotation file. SAF or GTF.
        annotation_format = 'GTF',
        # GTF-specific options.
        feature_type = 'exon',  # Feature type in GTF file.
        attribute_type = 'gene_id',  # Attribute type in GTF file.
        # Minimum number of overlapping bases in a read that is required for read assignment.
        # [default: 1]
        min_overlap = 1,
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/feature-counts'
