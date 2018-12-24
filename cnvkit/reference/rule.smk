# Please define a list of normal samples here.
NORMAL_SAMPLES = []

rule cnvkit_reference:
    input:
        target_coverage = expand('{sample}.targetcoverage.cnn', sample=NORMAL_SAMPLES),
        antitarget_coverage = expand('{sample}.antitargetcoverage.cnn', sample=NORMAL_SAMPLES),
        reference = 'reference.fasta',
    output:
        output_reference = 'reference.cnn',
    params:
        extra = '',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/reference'
