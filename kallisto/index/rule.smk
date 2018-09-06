rule kallisto_index:
    input:
        # Required input.
        'reference/Homo_sapiens.GRCh38.genome.fasta'
    output:
        # Required output.
        'reference/Homo_sapiens.GRCh38.genome.kallisto_idx'
    params:
        # Optional parameters. It omitted, default value will be used.
        extra = '',
    threads: 1
    log: 'logs/kallisto_index/Homo_sapiens.GRCh38.genome.log'
    wrapper:
        'http://dohlee-bio.info:9193/kallisto/index'
