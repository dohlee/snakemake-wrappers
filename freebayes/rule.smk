rule freebayes:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        bam = ['{sample}.bam'],
        reference = 'reference/Homo_sapiens_assembly38.fasta'
    output:
        # Required output.
        ['{sample}.vcf']
    params:
        extra = ''
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/freebayes'
