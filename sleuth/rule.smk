rule sleuth:
    input:
        expand(str(RESULT_DIR / '{sample}' / 'abundance.h5'), sample=SAMPLES)
    output:
        expand(str(RESULT_DIR / 'sleuth_deg.tsv'), sample=SAMPLES)
    params:
        # Manifest file that contains sample condition information.
        manifest = config['manifest'],
        # q-value cutoff to call differentially expressed genes.
        qval = 0.05,
    log: 'logs/sleuth.log'
    benchmark: 'benchmarks/sleuth.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sleuth/wrapper.R'
        
