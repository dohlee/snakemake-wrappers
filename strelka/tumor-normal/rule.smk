rule strelka_tumor_normal:
    input:
        # Required input.
        tumor = '{tumor_sample}.bam',
        normal = '{normal_sample}.bam',
        reference = 'reference/Homo_sapiens_assembly38.fasta'
    output:
        # Required output.
        snv = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.strelka.vcf.gz',
        indel = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.indels.strelka.vcf.gz'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        region = '',  # e.g. chr2:100-2000.
        call_regions = '',  # Here goes a predefined bed file.
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/strelka/tumor-normal'
