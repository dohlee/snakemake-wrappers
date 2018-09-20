rule somatic_sniper:
    input:
        # Required input.
        tumor_bam = '{tumor_sample}.bam',
        normal_bam = '{normal_sample}.bam',
        reference = 'reference/Homo_sapiens_assembly38.fasta'
    output:
        # Required output.
        'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.somatic_sniper.vcf',
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        mapping_quality_cutoff = 20,  # Recommended. (default 20)
        calling_quality_cutoff = 15,  # Recommended. (default 15)
        tumor_sample_name = '',
        normal_sample_name = ''
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/somatic-sniper'
