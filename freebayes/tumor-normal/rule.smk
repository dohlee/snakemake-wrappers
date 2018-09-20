rule freebayes_tumor_normal:
    input:
        # Required input.
        tumor_bam = '{tumor_sample}.sorted.bam',
        normal_bam = '{normal_sample}.sorted.bam',
        reference = 'reference/Homo_sapiens_assembly38.fasta',
        reference_index = 'reference/Homo_sapiens_assembly38.fasta.fai'
    output:
        # Required output.
        'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.freebayes.vcf'
    params:
        # Optional arguments. Omit if unneeded.
        extra = '',
        min_alternate_fraction = 0.03,  # Recommended.
        min_alternate_count = 2,  # Recommended.
        cnv_map = '',
    threads: 4
    log: 'logs/freebayes_tumor_normal/{tumor_sample}_vs_{normal_sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/freebayes/tumor-normal'
