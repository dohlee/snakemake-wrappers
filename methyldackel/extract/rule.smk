rule methyldackel_extract:
    input:
        bam = 'result/{sample}/{sample}.sorted.bam',
        bam_index = 'result/{sample}/{sample}.sorted.bam.bai'
        reference = 'reference/Homo_sapiens.GRCh38.genome.fasta'
    output:
        # NOTE: {sample}_CpG.meth.bedGraph will use --fraction option,
        # {sample}_CpG.counts.bedGraph will use --counts option, and
        # {sample}_CpG.logit.bedGraph will use --logit option.
        'result/{sample}/{sample}_CpG.bedGraph'
    threads: 4
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        min_depth = 10,
        mapping_quality_threshold = 10,
        sequencing_quality_threshold = 5,
    log: 'logs/methyldackel_extract/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/methyldackel/extract'
