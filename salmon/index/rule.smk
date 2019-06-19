rule salmon_index:
    input:
        'reference/Homo_sapiens.GRCh38.transcriptome.fasta'
    output:
        directory('reference/Homo_sapiens.GRCh38.transcriptome.salmon_idx')
    params:
        extra = '',
        # The size of k-mers that should be used for the quasi index.
        # Default: 31
        kmerLen = 31,
        # This flag will expect the input transcript fasta to be in 
        # GENCODE format, and will split the transcript name at the first '|' character.
        # These reduced names will be used in the output and when looking for these transcripts
        # in a gene to transcript GTF.
        # Default: False
        gencode = False,
        # This flag will disable the dfeault indexing behavior of discarding sequence-identical
        # duplicate transcripts. If this flag is passed, then duplicate transcripts that
        # appear in the input will be retained and quantified separately.
        # Default: False
        keepDuplicates = False,
        # [quasi index only] Build the index using a perfect hash rather than a dense hash.
        # This will require less memory (especially during quantification),
        # but will take longer to construct.
        # Default: False
        perfectHash = False,
        # Treat these sequences are decoys that may have sequence homologous to some known
        # transcript.
        # Default: False
        decoys = False,
        # The type of index to build; the only option is "quasi" in this version of salmon.
        # Default: quasi
        type_ = 'quasi',
    threads: 1
    log: 'logs/salmon/index/Homo_sapiens.GRCh38.transcriptome.log'
    benchmark: 'benchmarks/salmon/index/Homo_sapiens.GRCh38.transcriptome.log'
    wrapper:
        'http://dohlee-bio.info:9193/salmon/index'
