rule salmon_quant:
    input:
        # Required input.
        index = 'reference/Homo_sapiens.GRCh38.transcriptome.salmon_idx',
        fastq = ['data/{sample}.read1.fastq.gz', 'data/{sample}.read2.fastq.gz']
    output:
        # Required output.
        quant = 'result/{sample}/quant.sf',
        lib = 'result/{sample}/lib_format_counts.json',
    params:
        # Optional parameters. It omitted, default value will be used.
        library_type = 'A',
        extra = '',
    threads: 8
    log: 'logs/salmon_quant/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/salmon/quant'
