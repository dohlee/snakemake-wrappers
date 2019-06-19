rule bwameth:
    input:
        # Required input. Reference genome fasta file.
        reads = ['{sample}.fastq.gz'],
        reference = 'reference/hg38.fasta',
        bwameth_reference = 'reference/hg38.fasta.bwameth.c2t',
    output:
        # BAM output or SAM output is both allowed.
        # Note that BAM output will be automatically detected by its file extension,
        # and SAM output (which is bwa mem default) will be piped through `samtools view`
        # to convert SAM to BAM.
        'result/{sample}/{sample}.bam'
    params:
        extra = '',
        # Read-group to add to bam in same format as to bwa:
        # '@RG\tID:foo\tSM:bar'
        read_group = lambda wildcards: '@RG\tID:%s\tSM:%s\tPL:ILLUMINA' % (wildcards.sample, wildcards.sample),
        # flag alignments to this strand as not passing QC (0x200). Targetted BS-Seq libraries
        # are often to a single strand, so we can flag them as QC failures.
        # Note f == OT, r == OB. Likely, this will be 'f' as we will expect reads to align
        # to the original-bottom (OB) strand and will flag as failed those aligning to the
        # forward, or original top (OT).
        # Default: False
        set_as_failed = False
        # Fastq files have 4 lines of read1 followed by 4 lines of read 2.
        # e.g. seqtk mergepe output.
        # Default: False
        interleaved = False,
    threads: 8
    log: 'logs/bwameth/{sample}.log'
    benchmark: 'benchmarks/bwameth/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bwa-meth'
