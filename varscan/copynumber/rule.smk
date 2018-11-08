rule varscan_copynumber:
    input:
        # Required input.
        tumor_bam = '{tumor_sample}.sorted.bam',
        normal_bam = '{normal_sample}.sorted.bam',
        reference = 'reference/Homo_sapiens_assembly38.fasta'
    output:
        # Required output.
        # NOTE: Varscan can output variants in its own format such as output.snp, output.indel,
        # however, this wrapper forces VCF output for unity of variant calling pipelines.
        raw_copynumber_calls = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.copynumber'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Minimum base quality to count for coverage. (Default 20)
        min_base_qual = 20,
        # Minimum read mapping quality to count for coverage. (Default 20)
        min_map_qual = 20,
        # Minimum coverage threshold for copynumber segments. (Default 20)
        min_coverage = 20,
        # Minimum number of consecutive bases to report a segment. (Default 10)
        min_segment_size = 10,
        # Max size before a new segment is made. (Default 100)
        max_segment_size = 100,
        # P-value threshold for significant copynumber change-point. (Default 0.01)
        p_value = 0.01,
        # The normal/tumor input data ratio for copynumber adjustment. (Default 1.0)
        data_ratio = 1.0,
        pileup_quality_cutoff = 1,  # Recommended.  (Default 1)
        output_prefix = lambda wildcards, input: f'result/{wildcards.tumor_sample}/{wildcards.tumor_sample}_vs_{wildcards.normal_sample}',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/varscan/copynumber'
