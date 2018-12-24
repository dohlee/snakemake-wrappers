rule cnvkit_fix:
    input:
        target_coverage = '{tumor}.targetcoverage.cnn',
        antitarget_coverage = '{tumor}.antitargetcoverage.cnn',
        reference = 'reference.cnn',
    output:
        '{tumor}.cnr'
    params:
        extra = '',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/fix'
