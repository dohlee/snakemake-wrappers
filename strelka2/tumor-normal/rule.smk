rule strelka2_tumor_normal:
    input:
        # Required input.
        tumor = '{tumor_sample}.sorted.bam',
        tumor_index = '{tumor_sample}.sorted.bam.bai',
        normal = '{normal_sample}.sorted.bam',
        normal_index = '{normal_sample}.sorted.bam.bai',
        reference = 'reference/Homo_sapiens_assembly38.fasta',
        reference_index = 'reference/Homo_sapiens_assembly38.fasta.fai'
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
        'http://dohlee-bio.info:9193/strelka2/tumor-normal'
