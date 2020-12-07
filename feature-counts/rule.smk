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
        # Specify feature type(s) in a GTF annotation. If multiple types are provided,
        # they should be separated by ',' with no space in between.
        # 'exon' by default. Rows in the annotation with a matched feature will be
        # extracted and used for read mapping.
        feature_type = 'exon',
        # Specify attribute type in GTF annotation. 'gene_id' by default.
        # Meta-features used for read counting will be extracted from annotation
        # using the provided value.
        attribute_type = 'gene_id',
        # Extract extra attribute types from the provided GTF annotation and include
        # them in the counting output. These attribute types will not be used to group
        # features. If more than one attribute type is provided they should be
        # separated by comma.
        extra_attributes = False,
        # Provide a chromosome name alias file to match chr names in
        # annotation with those in the reads. This should be a two-column
        # comma-delimited text file. Its first column should include chr names
        # in the annotation and its second column should include chr names in the
        # reads. Chr names are case sensitive. No column header should be included
        # in the file.
        chromosome_name_alias = False,
        # Perform read counting at feature level (eg. counting reads for exons
        # rather than genes).
        feature_level_counting = False,
        # Assign reads to all their overlapping meta-features (or features if -f is
        # specified).
        multi_overlapping_reads = False,
        # Minimum number of overlapping bases in a read that is required for read
        # assignment. 1 by default. Number of overlapping bases is counted from both
        # reads if paired end. If a negative value is provided, then a gap of up to
        # specified size will be allowed between read and the feature that the read
        # is assigned to.
        min_overlap = 1,
        # Minimum fraction of overlapping bases in a read that is required for read
        # assignment. Value should be within range [0, 1]. 0 by default. Number of
        # overlapping bases is counted from both read if paired end. Both this option
        # and '--minOverlap' option need to be satisfied for read assignment.
        frac_overlap = 0.0,
        # Minimum fraction of overlapping bases in a feature that is required for
        # read assignment. Value should be with range [0, 1]. 0 by default.
        frac_overlap_feature = 0.0,
        # Assign reads to a meta-feature/feature that has the largest number of
        # overlapping bases.
        largest_overlap = False,
        # Maximum number of non-overlapping bases in a read (or a read pair)
        # that is allowed when being assigned to a feature. No limit is set by default.
        non_overlap = False,
        # Maximum number of non-overlapping bases in a feature that is allowed in 
        # read assignment. No limit is set by default.
        non_overlap_feature = False,
        # Reads are extended upstream by <int> bases from their 5' end.
        read_extension_5 = False,
        # Reads are extended upstream by <int> bases from their 3' end.
        read_extension_3 = False,
        # Reduce reads to their 5' most base or 3' most base. Read counting is then
        # performed based on the single base the read is reduced to.
        # e.g. --read2pos <5:3>
        read2pos = False,
        # Multi-mapping reads will also be counted. For a multi-mapping read,
        # all its reported alignments will be counted. The 'NH' tag in BAM/SAM
        # input is used to detect multi-mapping reads.
        count_multi_mapping = False,
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/feature-counts'
