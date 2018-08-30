rule bwa_mem:
    input:
        # Required input. Reference genome fasta file.
        reads = ['{sample}.read1.fastq.gz'],
        mates = ['{sample}.read2.fastq.gz'],
        # You may use any of {genome}.amb, {genome}.ann, {genome}.bwt,
        # {genome}.pac, {genome}.sa just to obligate snakemake to run `bwa index` first.
        reference = 'reference/hg38.bwt'
    output:
        # BAM output or SAM output is both allowed.
        # Note that BAM output will be automatically detected by its file extension,
        # and SAM output (which is bwa mem default) will be piped through `samtools view`
        # to convert SAM to BAM.
        'result/{sample}/{sample}.bam'
    params:
        # -M option marks secondary alignments.
        # You may need this if you use GATK downstream.
        extra = "-M " \
                # Read group annotation. Omit if unused.
                # NOTE: You should check the platform information of the read data!
                "-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
    threads: 8
    wrapper:
        'http://dohlee-bio.info:9193/bwa/mem'
