rule salmon_quant:
    input:
        # Required input.
        index = 'reference/Homo_sapiens.GRCh38.transcriptome.salmon_idx',
        reads = ['data/{sample}.read1.fastq.gz', 'data/{sample}.read2.fastq.gz']
    output:
        # Required output.
        quant = 'result/{sample}/quant.sf',
        lib = 'result/{sample}/lib_format_counts.json',
    params:
        extra = '',
        # Format string describing the library type.
        # Please refer to: https://salmon.readthedocs.io/en/latest/library_type.html
        # Default: 'A' (Automatically detect library type.)
        libType = 'A',
        # Perform sequence-specific bias correction.
        # Default: False
        seqBias = False,
        # Perform fragment GC bias correction.
        # Default: False
        gcBias = False,
        # This option sets the prior probability that an alignment that disagrees with
        # the specified library type (--libType) results from the true fragment origin.
        # Setting this to 0 specifies that alignments that disagree with the library type
        # are no less likely than those that do.
        # Default: False
        incompatPrior = False,
        # File containing a mapping of transcripts to genes. If this file is provided salmon
        # will output both quant.sf and quant.genes.sf files, where the latter contains 
        # aggregated gene-level abundance estimates. The transcript to gene mapping should be
        # provided as either a GTF file, or a in a simple tab-delimited format where each line
        # contains the name of a transcript and the gene to which it belongs separated by a tab.
        # The extension of the file is used to determine how the file should be parsed.
        # Files ending in '.gtf', '.gff' or '.gff3' are assumed to be in GTF format; files with
        # any other extension are assumed to be in the simple format. In GTF / GFF format, the 
        # "transcript_id" is assumed to contain the transcript identifier and the "gene_id" is
        # assumed to contain the correspondin gene identifier.
        # Default: False
        geneMap = False,
        # If you're using Salmon on a metagenomic dataset, consider setting this flag to disable
        # parts of the abundance estimation model that make less sense for metagenomic data.
        # Default: False
        meta = False,
        # NOTE: The rest of advanced options are omitted here.
        # so please refer to the documentation for whole options.
    threads: 8
    log: 'logs/salmon/quant/{sample}.log'
    benchmark: 'benchmarks/salmon/quant/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/salmon/quant'
