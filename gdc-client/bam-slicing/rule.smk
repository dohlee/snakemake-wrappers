rule bam_slicing:
    input:
        token = '',
    output:
        '{sample}.bam',
    params:
        # Required parameter.
        # GDC UUID for BAM file to slice.
        uuid = '',
        # Specify comma-separated region(s) if needed.
        region = 'chr2:100000-200000,chr3:100000-200000',
        # Specify comma-separated HGNC/GENCODE v22 gene name(s) if needed.
        gencode = 'BRCA1,DNMT3A',
    threads: 1
    log: 'logs/bam_slicing/{sample}.log'
    benchmark: 'benchmarks/bam_slicing/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/gdc-client/bam-slicing'
