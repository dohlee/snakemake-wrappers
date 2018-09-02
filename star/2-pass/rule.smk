rule star_2_pass:
    input:
        # Required input.
        # NOTE: Reference genome should be uncompressed.
        reads = ['{sample}.fastq.gz'],
        star_index = directory('reference/star_index')
    output:
        # There is no need to output sam or unsorted bam file!
        # So this wrapper includes '--outSAMtype BAM Unsorted' option by default.
        'result/{sample}/{sample}.sorted.bam'
    threads: 1
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = '',
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or alignReads step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        # NOTE: Make sure the annotation file is unzipped!
        sjdb_gtf_file = '',
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally = (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = 100,
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = '',
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = '',
    log: 'logs/star/star/2-pass/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/star/2-pass'
