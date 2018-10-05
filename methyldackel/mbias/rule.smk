rule methyldackel_mbias:
    input:
        bam = 'result/{sample}/{sample}.sorted.bam',
        bam_index = 'result/{sample}/{sample}.sorted.bam.bai',
        reference = 'reference/Homo_sapiens.GRCh38.genome.fasta',
    output:
        'result/{sample}/{sample}.mbias.tsv',
    threads: 4
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        mapping_quality_threshold = 10,
        sequencing_quality_threshold = 5,
    log: 'logs/methyldackel_mbias/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/methyldackel/mbias'
