rule star_genome_generate:
    input:
        # Required input.
        # NOTE: Reference genome should be uncompressed.
        reference = 'reference/Homo_sapiens.GRCh38.genome.fasta'
    output:
        index_directory = directory('reference/star_index/')
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
    log: 'logs/star/genome_generate/_.log'
    wrapper:
        'http://dohlee-bio.info:9193/star/genome-generate'
