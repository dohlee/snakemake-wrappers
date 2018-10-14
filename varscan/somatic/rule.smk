rule varscan_somatic:
    input:
        # Required input.
        tumor_bam = '{tumor_sample}.sorted.bam',
        normal_bam = '{normal_sample}.sorted.bam',
        reference = 'reference/Homo_sapiens_assembly38.fasta'
    output:
        # Required output.
        # NOTE: Varscan can output variants in its own format such as output.snp, output.indel,
        # however, this wrapper forces VCF output for unity of variant calling pipelines.
        snv = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.vcf',
        indel = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.indels.varscan.vcf'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Minimum frequency to call a heterozygote. (default 0.10)
        min_var_freq = 0.03,  # Recommended.
        # If set to 1, removes variants with >90% strand bias.
        strand_filter = 1,  # Recommended.
        pileup_quality_cutoff = 20,  # Recommended. (default 20)
        region = '',
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/varscan/somatic'
