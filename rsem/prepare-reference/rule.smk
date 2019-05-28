rule rsem_prepare_reference:
    input:
        fasta = 'YOUR_REFERENCE_FASTA'
        # You cannot provide both GTF and GFF3 files.
        gtf = 'ANNOTATION_GTF'
        gff3 = 'ANNOTATION_GFF3'
    output:
        'REFERENCE.transcripts.fa'
    params:
        extra = '',
        # Give a comma-separated list of transcript categories,
        # e.g. "mRNA,rRNA". Only transcripts that match the pattern
        # will be extracted.
        # Default: "mRNA"
        gff3_rna_patterns = 'mRNA',
        # Give a comma-separated list of trusted sources,
        # e.g. "ENSEMBL,HAVANA". Only transcripts coming from
        # these sources will be extracted. If this option is off,
        # all sources are accepted.
        # Default: False
        trusted_sources = False,
        # Use information from given file to map from transcript
        # (isoform) ids to gene ids. Each line of the file should be
        # of the form:
        #
        # gene_id transcript_id
        #
        # with the two fields separated by a tab character.
        # If you are using a GTF file for the "UCSC Genes" gene set from
        # the UCSC Genome Browser, then the "knownIsoforms.txt" file
        # (obtained from the "Downloads" section of the UCSC Genome Browser
        # site) is of this format.
        # If this option is off, then the mapping of isoforms to genes depends
        # on whether the '--gtf' option is specified. If '--gtf' is specified,
        # then RSEM uses the "gene_id" and "transcript_id" attributes in the GTF
        # file. Otherwise, RSEM assumes that each sequence in the reference
        # sequence files is a separate gene.
        # Default: False
        transcript_to_gene_map = False,
        # Use information from given file to provide gene_id and transcript_id
        # information for each allele-specific transcript. Each line of the file
        # should be of the form:
        #
        # gene_id transcript_id allele_id
        #
        # with the fields separated by a tab character.
        # This option is designed for quantifying allele-specific expression.
        # It is only valid if '--gtf' option is not specified. allele_id should be
        # the sequence names presented in the Multi-FASTA-formatted files.
        # Default: False
        allele_to_gene_map = False,
        # Supress the output of logging information.
        # Default: False
        quiet = False,
    threads: 4
    log: 'logs/rsem_prepare_reference/_.log'
    wrapper:
        'http://dohlee-bio.info:9193/rsem/prepare-reference'

