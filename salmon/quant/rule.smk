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
        # Discard orphan mappings in quasi-mapping mode. If this flag is passed then only paired
        # mappings will be considered toward quantification estimates. The default behavior is to
        # consider orphan mappings if no valid paired mappings exist. This flag is independent
        # of the option to write the orphaned mappings to file (--writeOrphanLinks).
        # Default: False
        discardOrphansQuasi = False,
        # Validate mappings using alignment-based verification. If this flag is passed, quasi-mappings
        # will be validated to ensure that they could give rise to a reasonable alignment before they are
        # further used for quantification.
        # Default: False
        validateMappings = False,
        # The amount of slack allowed in the quasi-mapping consensus mechanism. Normally, a transcript
        # must cover all hits to be considered for mapping. If this is set to a fraction, X, greater than
        # 0 (and in [0,1)), then a transcript can fail to cover up to (100 * X)% of the hits before it is
        # discounted as a mapping candidate. The default value of this option is 0.2 if --validateMappings
        # is given and 0 otherwise.
        # Default: False
        consensusSlack = False,
        # The fraction of the optimal possible alignment score that a mapping must achieve in order to be
        # considered "valid" --- should be in (0, 1].
        # Salmon Default 0.65 and Alevin Default 0.87
        # Default: False
        minScoreFraction = False,
        # Sets the maximum allowable MMP extension when collecting suffix array intervals to be used in
        # chaining. This prevents MMPs from becoming too long, and potentially masking intervals that would
        # uncover other goo quasi-mappings for the read. This heuristic mimics the idea of the
        # maximum mappable safe prefix (MMSP) in selective alignment. Setting a smaller value will potentially
        # allow for more sensitive, but slower, mapping.
        # Default: 7
        maxMMPExtension = 7,
        # The value given to a match between read and reference nucleotides in an alignment.
        # Default: 2
        ma = 2,
        # The value given to a mis-match between read and reference nucleotides in an alignment.
        # Default: -4
        mp = -4,
        # The value given to a gap opening in an alignment.
        # Default: 4
        go = 4,
        # The value given to a gap extension in an alignment.
        # Default: 2
        ge = 2,
        # The value used for the bandwidth passed to ksw2. A smaller bandwidth can make the alignment
        # verification run more quickly, but could possibly miss valid alignments.
        # Default: 15
        bandwidth = 15,
        # Allow dovetailing mappings.
        # Default: False
        allowDovetail = False,
        # Attempt to recover the mates of orphaned reads. This uses edlib for orphan recovery, and so
        # introduces some computational overhead, but it can improve sensitivity.
        # Default: False
        recoverOrphans = False,
        # Set flags to mimic parameters similar to Bowtie2 with --no-discordant and --no-mixed flags.
        # This increases disallows dovetailing reads, and discards orphans. Note, this does not impose
        # the very strict parameters assumed by RSEM+Bowtie2, like gapless alignments. For that behavior,
        # use the --mimicStrictBT2 flag below.
        # Default: False
        mimicBT2 = False,
        # Set flags to mimic the very strict parameters used by RSEM+Bowtie2. This increases --minScoreFraction
        # to 0.8, disallows dovetailing reads, discards orphans, and disallows gaps in alignments.
        # Default: False
        mimicStrictBT2 = False,
        # Instead of weighting mappings by their alignment score, this flag will discard any mappings with
        # sub-optimal alignment score. The default option of soft-filtering (i.e. weighting mappings by
        # their alignment score) usually yields slightly more accurate abundance estimates but this flag may be
        # desirable if you want more accurate 'naive' equivalence classes, rather than range factorizd equivalence
        # classes.
        hardFilter = False,
        # If this option is provided, then the quasi-mapping results will be written out in SAM-compatible
        # format. By default, output will be directed to stdout, but an alternative file name can be provided
        # instead.
        # Default: False (when specified, '-')
        writeMappings = False,
        # Force hits gathered during quasi-mapping to be "consistent" (i.e. co-linear and approximately the right
        # distance apart).
        # Default: False
        consistentHits = False,
        # Number of bootstrap samples to generate. Note: This is mutually exclusive with Gibbs sampling.
        # Default: False (0)
        numBootstraps = False,
        # NOTE: The rest of advanced options are omitted here
        # so please refer to the documentation for whole options.
    threads: 8
    log: 'logs/salmon/quant/{sample}.log'
    benchmark: 'benchmarks/salmon/quant/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/salmon/quant'
