# Please define a list of all samples here.
SAMPLES = []

rule cnvkit_autobin:
    input:
        bams = expand('{sample}.bam', sample=SAMPLES),
        targets = 'SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.bed',
        access = 'access-5k-mappable.hg19.bed',
    output:
        'SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.target.bed',
        'SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.antitarget.bed',
    params:
        extra = '',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/autobin'
