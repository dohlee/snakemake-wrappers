rule star_alignreads:
    input:
        # Required input.
        reads = ['{sample}.fastq.gz'],
        star_index = directory('reference/star_index')
    output:
        # There is no need to output sam or unsorted bam file!
        # So this wrapper includes '--outSAMtype BAM SortedByCoordinate' option by default.
        genome_alignment = 'result/{sample}/{sample}.sorted.bam',
        # Aligning to transcriptome is optional.
        # NOTE: You should set quantMode to 'TranscriptomeSAM' to get 
        # this output.
        transcriptome_alignment = 'result/{sample}/{sample}.bam'
    threads: 1
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = '',
        # Random number generator seed.
        # Default: 777
        run_rng_seed = 777,
        ### Quantification and Annotations
        # Types of quantification requested
        # False : none
        # TranscriptomeSAM: output SAM/BAM alignments to transcriptome into a separate file.
        # GeneCounts: count reads per gene
        # Default: False
        quant_mode = False,
        # Mode of shared memory usage for the genome files.
        # Only used with --runMode alignReads.
        # Default: NoSharedMemory
        genome_load = 'NoSharedMemory',
        # Path to the files with genomic coordinates
        # (chr <tab> start <tab> end <tab> strand)
        # for the splice junction introns. Multiple files can be supplied
        # and will be concatenated.
        # Default: False
        sjdb_file_chr_start_end = False,
        # It is recommended to give gene annotation file in at least one of
        # genomeGeneration or alignReads step.
        # Also it is recommened to use GENCODE annotation file, either of GTF of GFF3 file.
        # Don't forget to specify sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        # in case of GFF3 annotation file.
        # NOTE: Make sure the annotation file is unzipped!
        sjdb_gtf_file = '',
        # If genome is from UCSC, and annotation is from ENSEMBL,
        # use sjdb_gtf_chr_prefix = 'chr'
        sjdb_gtf_chr_prefix = '',
        # Length of the donor/acceptor sequence on each side of the junctions,
        # ideally, (maxReadLength - 1).
        # In most cases, the default value of 100 will work well.
        sjdb_overhang = 100,
        # If you use GFF3 annotation file,
        # use sjdb_gtf_tag_exon_parent_transcript = 'Parent'
        sjdb_gtf_tag_exon_parent_transcript = 'transcript_id',
        ### Output filtering
        # Type of filtering.
        # Normal: standard filtering using only current alignment.
        # BySJout: keep only those reads that contain junctions that passed 
        # filtering into SJ.out.tab
        # Default: Normal
        out_filter_type = 'Normal',
        # The score range below the maximum score for multimapping alignments.
        # Default: 1
        out_filter_multimap_score_range = 1,
        # Maximum number of loci the read is allowed to map to. Alignments (all of them)
        # will be output only if the read maps to no more loci than this value.
        # Otherwise no alignments will be output, and the read will be counted as
        # "mapped to too many loci" in the Log.final.out.
        # Default: 10
        out_filter_multimap_n_max = 10,
        # Alignment will be output only if it has no more mismatches than this value.
        # Default: 10
        out_filter_mismatch_n_max = 10,
        # Alignment will be output only if its ratio of mismatches to *mapped* length
        # is less than or equal to this value.
        # Default: 0.3
        out_filter_mismatch_n_over_l_max = 0.3,
        # Alignment will be output only if its ratio of mismatches to *read* length
        # is less than or equal to this value.
        # Default: 1.0
        out_filter_mismatch_n_over_read_l_max = 1.0,
        # Alignment will be output only if its score is higher than or equal
        # to this value.
        # Default: 0
        out_filter_score_min = 0,
        # Same as outFilterScoreMin, but normalized to read length.
        # (sum of mates' lengths for paired-end reads)
        # Default: 0.66
        out_filter_score_min_over_l_read = 0.66,
        # Alignment will be output only if the number of matched bases is higher than or equal to this value.
        # Default: 0
        out_filter_match_n_min = 0,
        # Same as outFilterMatchNmin, but normalized to the read length.
        # (sum of mates' lengths for pared-end reads)
        # Default: 0.66
        out_filter_match_n_min_over_l_read = 0.66,
        # Filter alignment using their motifs
        # None: no filtering
        # RemoveNoncanonical: filter out alignments that contain non-canonical junctions
        # RemoveNoncanonicalUnannotated: filter out alignments that contain non-canonical
        #   unannotated junctions when using annotated splice junctions database.
        #   The annotated non-canonical junctions will be kept.
        # Default: None
        out_filter_intron_motifs = 'None',
        # Filter alignments
        # RemoveInconsistentStrands: remove alignments that have junctions with
        # inconsistent strands. 
        # None: no filtering
        # Default: RemoveInconsistentStrands
        out_filter_intron_strands = 'RemoveInconsistentStrands',
        ### Output filtering: Splice junctions
        # Which reads to consider for collapsed splice junctions output.
        # All: all reads, unique- and multi-mappers
        # Unique: uniquely mapping reads only
        # Default: All
        out_sj_filter_reads = 'All',
        # Minimum overhang length for splice junctions on both sides for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and /GT/AT motif.
        # -1 means no output for that motif.
        # This does not apply to annotated junctions
        # Default: '30 12 12 12'
        out_sj_filter_overhang_min = '30 12 12 12',
        # Minimum uniquely mapping read count per junction for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and GT/AT motif.
        # -1 means no output for that motif.
        # Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
        # conditions are satisfied.
        # This does not apply to annotated junctions.
        # Default: '3 1 1 1'
        out_sj_filter_count_unique_min = '3 1 1 1',
        # Minimum total (multi-mapping+unique) read count per junction for:
        # (1) non-canonical motifs,
        # (2) GT/AG and CT/AC motif,
        # (3) GC/AG and CT/GC motif,
        # (4) AT/AC and GT/AT motif.
        # -1 means no output for that motif.
        # Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin
        # conditions are satisfied.
        # This does not apply to annotated junctions.
        # Default: '3 1 1 1'
        out_sj_filter_count_total_min = '3 1 1 1',
        # Minimum allowed distance to other junctions' donor/acceptor.
        # This does not apply to annotated junctions.
        # Default: '10 0 5 10'
        out_sj_filter_dist_to_other_sj_min = '10 0 5 10',
        # Maximum gap allowed for junctions supported by 1,2,3,,,N reads.
        # i.e. by default junctions supported
        # by 1 read can have gaps <= 50000,
        # by 2 reads: <= 100000,
        # by 3 reads: <= 200000,
        # by >=4 reads any gap <= alignIntronMax
        # This does not apply to annotated junctions.
        # Default: '50000 100000 200000'
        out_sj_filter_intron_max_vs_read_n = '50000 100000 200000',
        ### Alignments and seeding
        # Minimum intron size: genomic gap is considered intron if
        # its length>=alignIntronMin, otherwise it is considered deletion.
        # Default: 21
        align_intron_min = 21,
        # maximum intron size: if 0, max intron size will be determined by
        # (2^winBinNbits)*winAnchorDistNbins
        # Default: 0
        align_intron_max = 0,
        # Minimum overhang (i.e. block size) for spliced alignments.
        # Default: 5
        align_sj_overhang_min = 5,
    log: 'logs/star_alignreads/{sample}.log'
    benchmark: 'benchmarks/star_alignreads/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/star/alignreads'
