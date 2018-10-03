rule varscan_process_somatic:
    input:
        # Required input.
        'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.vcf',
        # OR
        #'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.indels.varscan.vcf',
    output:
        # Required output.
        germline_hc = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.Germline.hc.vcf',
        germline = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.Germline.vcf',
        loh_hc = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.LOH.hc.vcf',
        loh = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.LOH.vcf',
        somatic_hc = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.Somatic.hc.vcf',
        somatic = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.varscan.Somatic.vcf',
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Minimum variant allele frequency in tumor. (default 0.10)
        min_tumor_freq = 0.10,
        # Maximum variant allele frequency in normal. (default 0.05)
        max_normal_freq = 1,
        # P-value for high-confidence calling. (default 0.07)
        p_value = 0.07,
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/varscan/process-somatic'
