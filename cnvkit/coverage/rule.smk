rule cnvkit_coverage:
    input:
        bam = '{sample}.bam',
        targets = 'SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.target.bed',
        antitargets = 'SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.antitarget.bed',
    output:
        target_coverage = '{sample}.targetcoverage.cnn',
        antitarget_coverage = '{sample}.antitargetcoverage.cnn',
    params:
        extra = '',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/coverage'
