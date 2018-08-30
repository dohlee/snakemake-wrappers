rule bwa_meth:
    input:
        # Required input. Reference genome fasta file.
        reads = ['{sample}.read1.fastq.gz'],
        mates = ['{sample}.read2.fastq.gz'],
        reference = 'reference/hg38.fasta'
    output:
        # BAM output or SAM output is both allowed.
        # Note that BAM output will be automatically detected by its file extension,
        # and SAM output (which is bwa mem default) will be piped through `samtools view`
        # to convert SAM to BAM.
        'result/{sample}/{sample}.bam'
    params:
        # Read group annotation. Omit if unused.
        # NOTE: You should check the platform information of the read data!
        extra = "--read-group '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",

    threads: 8
    wrapper:
        'http://dohlee-bio.info:9193/bwa/meth'
