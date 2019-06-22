from pathlib import Path

rule wasabi:
    input:
        # For example,
        expand(str(RESULT_DIR / '{sample}' / 'quant.sf'), sample=config['samples']),
    output:
        expand(str(RESULT_DIR / '{sample}' / 'abundance.h5'), sample=config['samples']),
    log: 'logs/wasabi/{sample}.log'
    benchmark: 'benchmarks/wasabi/{sample}.benchmark'
    script:
        'http://dohlee-bio.info:9193/wasabi/wrapper.R'

