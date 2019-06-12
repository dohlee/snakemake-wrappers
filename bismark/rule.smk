rule bismark_se:
    input:
        fastq = 'data/{sample}.trimmed.fastq.gz',
        reference_dir = directory('reference'),
        bisulfite_genome_dir = directory('reference/Bisulfite_Genome'),
    output:
        'result/se/{sample}/{sample}.bismark.bam',
        'result/se/{sample}/{sample}.bismark_report.txt',
    params:
        extra = '',
        # Sets single-end mapping mode explicitly giving a list of file names as <list>.
        # The filenames may be provided as a comma [,] or colon [:] separated list.
        # Default: False
        single_end = False,
        # The query input files (specified as <mate1>,<mate2> or <singles>) are FASTQ files
        # (usually having extension .fg or .fastq). This is the default. See also --solexa-quals.
        # Default: True
        fastq = True,
        # The query input files (specified as <mate1>,<mate2> or <singles>) are FASTA files
        # (usually having extensions .fa, .mfa, .fna or similar). All quality values are
        # assumed to be 40 on the Phred scale. FASTA files are expected to contain both the
        # read name and the sequence on a single line (and not spread over several lines).
        # Default: False
        fasta = False,
        # Skip (i.e. do not align) the first <int> reads or read pairs from the input.
        # Default: False
        skip = False,
        # Only aligns the first <int> reads or read pairs from the input.
        # Default: False (no limit)
        upto = False,
        # FASTQ qualities are ASCII chars equal to the Phred quality plus 33.
        # Default: True
        phred33_quals = True,
        # FASTQ qualities are ASCII chars equal to the Phred quality plus 64.
        # Default: False
        phred64_quals = False,
        # The full path to the Bowtie 2 installation on your system.
        # If not specified it is assumed that Bowtie 2 is in the PATH.
        # Default: False
        path_to_bowtie2 = False,
        # The full path to the HISAT2 installation on your system.
        # If not specified it is assumed that HISAT2 is in the PATH.
        path_to_hisat2 = False,
        #
        # Alignment
        #
        # Sets the number of mismatches to allowed in a seed alignment during multiseed alignment.
        # Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower)
        # but increases sensitivity.
        # Default 0. This option is only available for Bowtie 2 (for Bowtie 1 see -n)
        N = 0,
        # Sets the length of the seed substrings to align during multiseed alignment.
        # Smaller values make aligment slower but more sensitive.
        # The --sensitive preset of Bowtie 2 is used by default, which sets -L to 20.
        # Maximum of L can be set to 32. THe length of the seed would effect the alignment speed
        # dramatically while the larger L, the faster the alignment.
        # This option is only available for Bowtie 2 (for Bowtie 1 see -l).
        # Default: False (--sensitive preset of Bowtie 2)
        L = False,
        # When calculating a mismatch penalty, always consider the quality value at the mismatched
        # position to be the highest possible, regardless of the actual value. i.e. input is treated
        # as though all quality values are high. This is also the default behavior when the input
        # doesn't specify quality values (e.g. in -f mode). This option is invariable and on by default.
        # Default: False,
        ignore_quals = False,
        # The minimum insert size for valid paired-end alignments. e.g. if -I 60 is specified and
        # a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with
        # a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied).
        # A 19-bp gap would not be valid in that case.
        # Default: 0
        minins = 0,
        # The maixmum insert size for valid paired-end alignments. e.g. if -X 100 is specified and
        # a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 
        # 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied)
        # A 61-bp gap would not be valid in that case.
        # Default: 500
        maxins = 500,
        # In this mode, it is not required that the entire read aligns from one end to the other.
        # Rather, some characters may be omitted ("soft-clipped") from the ends in order to achieve the 
        # greatest alignment score. For Bowtie 2, the match bonus --ma (default: 2) is used in this mode,
        # and the best possible alignment score is equal to the match bonus (--ma) times the length of the
        # read. This is mutually exclusive with end-to-end alignments. For HISAT2, it is currently not
        # exactly known how the best alignment is calculate.
        # Default: False
        local = False,
        # The sequencing library was constructed in a non strand-specific manner, alignments to all four
        # bisulfite strands will be reported. The current Illumina protocol for BS-Seq is directional, in
        # which case the strands complementary to the original strands are merely theoretical and should
        # not exist in reality. Specifying directional alignments (which is the default) will only run
        # 2 alignment threads to the original top (OT) or bottom (OB) strands in parallel and report these
        # alignments. This is the recommended option for strand-specific libraries.
        # Deafult: False
        non_directional = False,
        # This option maay be used for PBAT-Seq libraries (Post-Bisulfite Adapter Tagging; Kobayashi et al.,
        # PLoS Genetics, 2012). This is essentially the exact opposite of alignments in 'directional' mode,
        # as it will only launch two alignment threads to the CTOT and CTOB strands instead of the normal OT
        # and OB ones. Use the option only if you are ceratin that your libraries were constructed following
        # a PBAT protocol (if you don't know what PBAT-Seq is you should not specify this option).
        # The option --pbat works only for FastQ files (in both Bowtie and Bowtie 2mode) and using uncompressed
        # temporary files only).
        # Default: False
        pbat = False,
        # Suppress SAM header lines (starting with @). THis might be useful when very large input files are
        # split up into several smaller files to run concurrently and the output files are to be merged.
        # Default: False
        sam_no_hd = False,
        # Write out a Read Group tag to the resulting SAM/BAM file. This will write the following line
        # to the SAM header: @RG PL: ILLUMINA ID:SAMPLE SM:SAMPLE ; to set ID and SM see --rg_id and --rg_sample.
        # In addition each read receives an RG:Z:RG-ID tag.
        rg_tag = False,
        # Sets the ID field in the @RG header line.
        # The default is 'SAMPLE'
        # Default: False (SAMPLE)
        rg_id = False,
        # Sets the SM field in the @RG header line; can't be wet without setting --rg_id as well.
        # The default is 'SAMPLE'
        # Default: False (SAMPLE)
        rg_sample = False,
        # Write all reads that could not be aligned to a file in the output directory. Written reads
        # will appear as they did in the input, without any translation of quality values that may have
        # taken place within Bowtie or Bismark. Paired-end reads will be written to two parallel files
        # with _1 and _2 inserted in their filenames, i.e. _unmapped_reads_.txt and unmapped_reads_2.txt.
        # Reads with more than one valid alignment with the same number of lowest mismatches
        # (ambiguous mapping) are also written to _unmapped_reads.txt unless the option --ambiguous is
        # specified as well.
        # Default: False
        unmapped = False,
        # Write all reads which produce more than one valid alignment with the same number of lowest mismatches
        # or other reads that fail to align uniquely to a file in the output directory.
        # Written reads will appear as they did in the input, without any of the translation of quality
        # values that may have taken place within Bowtie or Bismark. Paired-end reads will be written
        # parallel files with _1 and _2 inserted in their filenames, i.e. _ambiguous_reads_1.txt and
        # _ambiguous_reads_2.txt. These reads are not written to the file specified with --un.
        # Default: False
        ambiguous = False,
        # Write temporary files to this directory instead of into the same directory as the input files. If
        # the specified folder does not exist, Bismark will attempt to create it first. The path to the
        # temporary folder can be either relative or absolute.
        # Default: False
        temp_dir = False,
        # Optionally outputs an extra column specifying the number of non-bisulfite mismatches a read
        # during the alignment step. This option is only available for SAM format. In Bowtie 2 context, this
        # value is just the number of actual non-bisulfite mismatches and ignores potential insertions
        # or deletions. The format for single-end reads and read 1 of paired-end reads is 'XA:Z:number of
        # mismatches' and 'XB:Z:number of mismathces' for read 2 of paired-end reads.
        # Default: False
        non_bs_mm = False,
        # Temporary bisulfite conversion files will be written out in a GZIP compressed form to save disk space.
        # This option is available for most alignment modes but is not available for paired-end FastA files.
        # This option might be somewhat slower than writin gout uncompressed files, but this awaits further
        # testing.
        # Default: False
        gzip = False,
        # The output will be written out in SAM format instead of the default BAM format. BIsmark will attempt
        # to use the path to samtools tha twas specified with '--samtools_path', or, if it hasn't been specified,
        # attempt to find Samtools in the PATH. If no installation of Samtools can be found, the SAM output
        # will be compressed with GZIP instead (yielding a .sam.gz output file).
        # Default: False
        sam = False,
        # Writes the output to a CRAM file instread of BAM. This requires the use of samtools 1.2 or higher.
        # Default: False
        cram = False,
        # CRAM output requires you to specify a reference genome as a single FastA file. If this single-FastA
        # reference file is not supplied explicitly it will be regenerated from the genome .fa sequence(s)
        # used for the Bismark run and written to a file called 'Bismark-genome_CRAM_reference.mfa' into the
        # output directory.
        # Default: False
        cram_ref = False,
        # The path to your samtools installation, e.g., /home/user/samtools/. Does not need to be specified
        # if samtools is in the PATH already.
        # Default: False
        samtools_path = False,
        # Prefixes <prefix> to the output filenames. Trailing dots will be replaced by a single one. For example,
        # '--prefix test' with 'file.fa' would result in the output file 'test.file.fq_bismark.sam' etc.
        # Default: False
        prefix = False,
        # Write all output to files starting with this base file name. For example, '--basename foo'
        # would result in the files 'foo.bam' and 'foo_SE_report.txt' (or its paired-end equivalent).
        # Takes precedence over --prefix. Be advised that you should not use this option in conjunction with
        # supplying list of files to be processed consecutively, as all output files will constantly overwrite
        # each other.
        # Default: False
        basename = False,
        # Only in paired-end SAM mode, uses the FLAG values used by Bismark v0.8.2 and before. In addition, this
        # options appends /1 and /2 to the read IDs for reads 1 and 2 relative to the input file.
        # Since both the appended read IDs and custom FLAG values may cause problems with some downstream tools
        # such as Picard, new defaults were implemented as of version 0.8.3.
        # Default: False
        old_flag = False,
        # For reads that have multiple alignments a random alignment is written out to a special file
        # ending in '.ambiguous.bam'. The alignments are in Bowtie2 format and do not any contain Bismark
        # specific entries such as the methylation call etc. These ambiguous BAM files are intended to be used
        # as coverage estimators for variant callers.
        # Default: False
        ambig_bam = False,
        # Calculates the mono- and di-nucleotide sequence composition of covered positions in the analysed
        # BAM file and compares it to the genomic average composition once alignments are complete by calling
        # 'bam2nuc'. Since this calculation may take a while, bam2nuc attempts to write the genomic sequence
        # composition into a file called 'genomic_nucleotide_frequencies.txt' inside the reference genome folder
        # so it can be re-used the next time round instead of calculating it once again. If a file
        # 'nucleotide_stats.txt' is found with the Bismark reports it will be automatically detected and used
        # for the Bismark HTML report. This option works only for BAM or CRAM files.
        # Default: False
        nucleotide_coverage = False,
        # This option will truncate read IDs at the first space or tab it encounters, which are sometimes used
        # to add comments to a FastQ entry (instread of replacing them with underscores (_) as is the default
        # behaviour). The option is deliberately somewhat cryptic ("I couldn't possibly comment"), as it
        # only becomes relevant when R1 and R2 of read pairs are mapped separately in single-end mode, and
        # then re-paired afterwards (the SAM format dictates that R1 and R2 have the same read ID). Paired-end
        # mapping already creates BAM files with identical read IDs. Formore information please see here:
        # https://github.com/FelixKruger/Bismark/issues/236.
        # Default: False
        icpc = False,
    threads: 6
    log: 'logs/bismark_se/{sample}.log'
    benchmark: repeat('benchmarks/bismark_se/{sample}.benchmark', 1)
    wrapper:
        'http://dohlee-bio.info:9193/bismark'

rule bismark_pe:
    input:
        fastq = ['data/{sample}.read1.trimmed.fastq.gz', 'data/{sample}.read2.trimmed.fastq.gz'],
        reference_dir = directory('reference'),
        bisulfite_genome_dir = directory('reference/Bisulfite_Genome'),
    output:
        'result/pe/{sample}/{sample}.bismark.bam',
        'result/pe/{sample}/{sample}.bismark_report.txt',
    params:
        extra = '',
        # Sets single-end mapping mode explicitly giving a list of file names as <list>.
        # The filenames may be provided as a comma [,] or colon [:] separated list.
        # Default: False
        single_end = False,
        # The query input files (specified as <mate1>,<mate2> or <singles>) are FASTQ files
        # (usually having extension .fg or .fastq). This is the default. See also --solexa-quals.
        # Default: True
        fastq = True,
        # The query input files (specified as <mate1>,<mate2> or <singles>) are FASTA files
        # (usually having extensions .fa, .mfa, .fna or similar). All quality values are
        # assumed to be 40 on the Phred scale. FASTA files are expected to contain both the
        # read name and the sequence on a single line (and not spread over several lines).
        # Default: False
        fasta = False,
        # Skip (i.e. do not align) the first <int> reads or read pairs from the input.
        # Default: False
        skip = False,
        # Only aligns the first <int> reads or read pairs from the input.
        # Default: False (no limit)
        upto = False,
        # FASTQ qualities are ASCII chars equal to the Phred quality plus 33.
        # Default: True
        phred33_quals = True,
        # FASTQ qualities are ASCII chars equal to the Phred quality plus 64.
        # Default: False
        phred64_quals = False,
        # The full path to the Bowtie 2 installation on your system.
        # If not specified it is assumed that Bowtie 2 is in the PATH.
        # Default: False
        path_to_bowtie2 = False,
        # The full path to the HISAT2 installation on your system.
        # If not specified it is assumed that HISAT2 is in the PATH.
        path_to_hisat2 = False,
        #
        # Alignment
        #
        # Sets the number of mismatches to allowed in a seed alignment during multiseed alignment.
        # Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower)
        # but increases sensitivity.
        # Default 0. This option is only available for Bowtie 2 (for Bowtie 1 see -n)
        N = 0,
        # Sets the length of the seed substrings to align during multiseed alignment.
        # Smaller values make aligment slower but more sensitive.
        # The --sensitive preset of Bowtie 2 is used by default, which sets -L to 20.
        # Maximum of L can be set to 32. THe length of the seed would effect the alignment speed
        # dramatically while the larger L, the faster the alignment.
        # This option is only available for Bowtie 2 (for Bowtie 1 see -l).
        # Default: False (--sensitive preset of Bowtie 2)
        L = False,
        # When calculating a mismatch penalty, always consider the quality value at the mismatched
        # position to be the highest possible, regardless of the actual value. i.e. input is treated
        # as though all quality values are high. This is also the default behavior when the input
        # doesn't specify quality values (e.g. in -f mode). This option is invariable and on by default.
        # Default: False,
        ignore_quals = False,
        # The minimum insert size for valid paired-end alignments. e.g. if -I 60 is specified and
        # a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with
        # a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied).
        # A 19-bp gap would not be valid in that case.
        # Default: 0
        minins = 0,
        # The maixmum insert size for valid paired-end alignments. e.g. if -X 100 is specified and
        # a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 
        # 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied)
        # A 61-bp gap would not be valid in that case.
        # Default: 500
        maxins = 500,
        # In this mode, it is not required that the entire read aligns from one end to the other.
        # Rather, some characters may be omitted ("soft-clipped") from the ends in order to achieve the 
        # greatest alignment score. For Bowtie 2, the match bonus --ma (default: 2) is used in this mode,
        # and the best possible alignment score is equal to the match bonus (--ma) times the length of the
        # read. This is mutually exclusive with end-to-end alignments. For HISAT2, it is currently not
        # exactly known how the best alignment is calculate.
        # Default: False
        local = False,
        # The sequencing library was constructed in a non strand-specific manner, alignments to all four
        # bisulfite strands will be reported. The current Illumina protocol for BS-Seq is directional, in
        # which case the strands complementary to the original strands are merely theoretical and should
        # not exist in reality. Specifying directional alignments (which is the default) will only run
        # 2 alignment threads to the original top (OT) or bottom (OB) strands in parallel and report these
        # alignments. This is the recommended option for strand-specific libraries.
        # Deafult: False
        non_directional = False,
        # This option maay be used for PBAT-Seq libraries (Post-Bisulfite Adapter Tagging; Kobayashi et al.,
        # PLoS Genetics, 2012). This is essentially the exact opposite of alignments in 'directional' mode,
        # as it will only launch two alignment threads to the CTOT and CTOB strands instead of the normal OT
        # and OB ones. Use the option only if you are ceratin that your libraries were constructed following
        # a PBAT protocol (if you don't know what PBAT-Seq is you should not specify this option).
        # The option --pbat works only for FastQ files (in both Bowtie and Bowtie 2mode) and using uncompressed
        # temporary files only).
        # Default: False
        pbat = False,
        # Suppress SAM header lines (starting with @). THis might be useful when very large input files are
        # split up into several smaller files to run concurrently and the output files are to be merged.
        # Default: False
        sam_no_hd = False,
        # Write out a Read Group tag to the resulting SAM/BAM file. This will write the following line
        # to the SAM header: @RG PL: ILLUMINA ID:SAMPLE SM:SAMPLE ; to set ID and SM see --rg_id and --rg_sample.
        # In addition each read receives an RG:Z:RG-ID tag.
        rg_tag = False,
        # Sets the ID field in the @RG header line.
        # The default is 'SAMPLE'
        # Default: False (SAMPLE)
        rg_id = False,
        # Sets the SM field in the @RG header line; can't be wet without setting --rg_id as well.
        # The default is 'SAMPLE'
        # Default: False (SAMPLE)
        rg_sample = False,
        # Write all reads that could not be aligned to a file in the output directory. Written reads
        # will appear as they did in the input, without any translation of quality values that may have
        # taken place within Bowtie or Bismark. Paired-end reads will be written to two parallel files
        # with _1 and _2 inserted in their filenames, i.e. _unmapped_reads_.txt and unmapped_reads_2.txt.
        # Reads with more than one valid alignment with the same number of lowest mismatches
        # (ambiguous mapping) are also written to _unmapped_reads.txt unless the option --ambiguous is
        # specified as well.
        # Default: False
        unmapped = False,
        # Write all reads which produce more than one valid alignment with the same number of lowest mismatches
        # or other reads that fail to align uniquely to a file in the output directory.
        # Written reads will appear as they did in the input, without any of the translation of quality
        # values that may have taken place within Bowtie or Bismark. Paired-end reads will be written
        # parallel files with _1 and _2 inserted in their filenames, i.e. _ambiguous_reads_1.txt and
        # _ambiguous_reads_2.txt. These reads are not written to the file specified with --un.
        # Default: False
        ambiguous = False,
        # Write temporary files to this directory instead of into the same directory as the input files. If
        # the specified folder does not exist, Bismark will attempt to create it first. The path to the
        # temporary folder can be either relative or absolute.
        # Default: False
        temp_dir = False,
        # Optionally outputs an extra column specifying the number of non-bisulfite mismatches a read
        # during the alignment step. This option is only available for SAM format. In Bowtie 2 context, this
        # value is just the number of actual non-bisulfite mismatches and ignores potential insertions
        # or deletions. The format for single-end reads and read 1 of paired-end reads is 'XA:Z:number of
        # mismatches' and 'XB:Z:number of mismathces' for read 2 of paired-end reads.
        # Default: False
        non_bs_mm = False,
        # Temporary bisulfite conversion files will be written out in a GZIP compressed form to save disk space.
        # This option is available for most alignment modes but is not available for paired-end FastA files.
        # This option might be somewhat slower than writin gout uncompressed files, but this awaits further
        # testing.
        # Default: False
        gzip = False,
        # The output will be written out in SAM format instead of the default BAM format. BIsmark will attempt
        # to use the path to samtools tha twas specified with '--samtools_path', or, if it hasn't been specified,
        # attempt to find Samtools in the PATH. If no installation of Samtools can be found, the SAM output
        # will be compressed with GZIP instead (yielding a .sam.gz output file).
        # Default: False
        sam = False,
        # Writes the output to a CRAM file instread of BAM. This requires the use of samtools 1.2 or higher.
        # Default: False
        cram = False,
        # CRAM output requires you to specify a reference genome as a single FastA file. If this single-FastA
        # reference file is not supplied explicitly it will be regenerated from the genome .fa sequence(s)
        # used for the Bismark run and written to a file called 'Bismark-genome_CRAM_reference.mfa' into the
        # output directory.
        # Default: False
        cram_ref = False,
        # The path to your samtools installation, e.g., /home/user/samtools/. Does not need to be specified
        # if samtools is in the PATH already.
        # Default: False
        samtools_path = False,
        # Prefixes <prefix> to the output filenames. Trailing dots will be replaced by a single one. For example,
        # '--prefix test' with 'file.fa' would result in the output file 'test.file.fq_bismark.sam' etc.
        # Default: False
        prefix = False,
        # Write all output to files starting with this base file name. For example, '--basename foo'
        # would result in the files 'foo.bam' and 'foo_SE_report.txt' (or its paired-end equivalent).
        # Takes precedence over --prefix. Be advised that you should not use this option in conjunction with
        # supplying list of files to be processed consecutively, as all output files will constantly overwrite
        # each other.
        # Default: False
        basename = False,
        # Only in paired-end SAM mode, uses the FLAG values used by Bismark v0.8.2 and before. In addition, this
        # options appends /1 and /2 to the read IDs for reads 1 and 2 relative to the input file.
        # Since both the appended read IDs and custom FLAG values may cause problems with some downstream tools
        # such as Picard, new defaults were implemented as of version 0.8.3.
        # Default: False
        old_flag = False,
        # For reads that have multiple alignments a random alignment is written out to a special file
        # ending in '.ambiguous.bam'. The alignments are in Bowtie2 format and do not any contain Bismark
        # specific entries such as the methylation call etc. These ambiguous BAM files are intended to be used
        # as coverage estimators for variant callers.
        # Default: False
        ambig_bam = False,
        # Calculates the mono- and di-nucleotide sequence composition of covered positions in the analysed
        # BAM file and compares it to the genomic average composition once alignments are complete by calling
        # 'bam2nuc'. Since this calculation may take a while, bam2nuc attempts to write the genomic sequence
        # composition into a file called 'genomic_nucleotide_frequencies.txt' inside the reference genome folder
        # so it can be re-used the next time round instead of calculating it once again. If a file
        # 'nucleotide_stats.txt' is found with the Bismark reports it will be automatically detected and used
        # for the Bismark HTML report. This option works only for BAM or CRAM files.
        # Default: False
        nucleotide_coverage = False,
        # This option will truncate read IDs at the first space or tab it encounters, which are sometimes used
        # to add comments to a FastQ entry (instread of replacing them with underscores (_) as is the default
        # behaviour). The option is deliberately somewhat cryptic ("I couldn't possibly comment"), as it
        # only becomes relevant when R1 and R2 of read pairs are mapped separately in single-end mode, and
        # then re-paired afterwards (the SAM format dictates that R1 and R2 have the same read ID). Paired-end
        # mapping already creates BAM files with identical read IDs. Formore information please see here:
        # https://github.com/FelixKruger/Bismark/issues/236.
        # Default: False
        icpc = False,
    threads: 6
    log: 'logs/bismark_pe/{sample}.log'
    benchmark: repeat('benchmarks/bismark_pe/{sample}.benchmark', 1)
    wrapper:
        'http://dohlee-bio.info:9193/bismark'
