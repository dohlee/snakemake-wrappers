rule salmon_index:
    input:
        'reference/Homo_sapiens.GRCh38.transcriptome.fasta'
    output:
        directory('reference/Homo_sapiens.GRCh38.transcriptome.salmon_idx')
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
    threads: 1
    log: 'logs/salmon/index/Homo_sapiens.GRCh38.transcriptome.log'
    wrapper:
        'http://dohlee-bio.info:9193/salmon/index'
