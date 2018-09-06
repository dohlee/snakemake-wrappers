rule cufflinks:
    input:
        # Required input.
        alignment = ['{sample}.bam'],
        # Optional input. Omit if unneeded.
        annotation = 'annotation/gencode.v28.gtf.gz'
    output:
        # Required output.
        'result/{sample}/transcripts.gtf',
        'result/{sample}/genes.fpkm_tracking',
        'result/{sample}/isoforms.fpkm_tracking'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        random_seed = 1,

    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/cufflinks'
