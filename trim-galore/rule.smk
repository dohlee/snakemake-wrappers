rule trim_galore_se:
    input:
        '{sample}.fastq.gz'
    output:
        '{sample}.trimmed.fastq.gz'
    params:
        extra = '',
        # Trim low-quality ends from reads in addition to adapter removal. For RRBS samples,
        # quality trimming will be performed first, and adapter trimming is carried in a second round.
        # Other files are quality and adapter trimmed in a single pass. The algorithm is the same
        # as the one used by BWA (Subtract INT from all qualities; compute partial sums from all
        # indices to the end of the sequence; cut sequence at the index at which the sum is minimal).
        # Default: 20
        quality = 20,
        # Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding)
        # for quality trimming.
        # Default: True
        phred33 = True,
        # Instructs Cutadapt to use ASCII+64 quality scores as Phred scores (Illumina 1.5 encoding) for
        # quality trimming.
        # Default: False
        phred64 = False,
        # Run FastQC in the default mode on the FastQ file once trimming is complete.
        # Default: False
        fastqc = False,
        # Passes extra arguments to FastQC. If more than one argument is to be passed to FastQC they must
        # be in the form "arg1 arg2 etc.". An example would be: --fastqc_args "--nogroup --outdir /home/".
        # Passing extra arguments will automatically invoke FastQC, so --fastqc does not have to be specified
        # separately.
        # Default: False
        fastqc_args = False,
        # Adapter sequence to be trimmed. If no specified explicitly, Trim Galore will try to auto-detect
        # whether the Illumina universal, Nextera transposase or Illumina small RNA adapter sequence was used.
        # Also see '--illumina', '--nextera' and '--small_rna'. If no adapter can be detected within the
        # first 1 million sequences of the first file specified Trim Galore defaults to '--illumina'.
        # A single base may also be given as e.g. -a A{10}, to be expanded to -a AAAAAAAAAA.
        # Default: False
        adapter = False,
        # Optional adapter sequence to be trimmed off read 2 of paired-end files. This option requires
        # '--paired' to be specified as well. If the libraries to be trimmed are smallRNA then a2 will
        # be set to the Illumina small RNA 5' adapter automatically (GATCGTCGGACT). A single base may
        # also be given as e.g. -a2 A{10}, to be expanded to -a2 AAAAAAAAAA.
        # Default: False
        adapter2 = False,
        # Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter 'AGATCGGGAGC'
        # instead of the default auto-detection of adapter sequence.
        # Default: False
        illumina = False,
        # Adapter sequence to be trimmed is the first 12bp of the Nextera adapter 'CTGTCTCTTATA' instead of
        # default auto-detection of adapter sequence.
        # Default: False
        nextera = False,
        # Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter
        # 'TGGTTCTCGG' instead of the default auto-detection of adapter sequence. Selecting to trim
        # smallRNA adapters will also lower the --length value to 18bp. If the smallRNA libraries
        # are paired-end then a2 will be set to the Illumina small RNA 5' adapter automatically
        # (GATCGTCGGACT) unless -a 2 had been defined explicitly.
        # Default: False
        small_rna = False,
        # Discard reads that are longer than <INT> bp after trimming. This is only advised for smallRNA
        # sequencing to remove non-small RNA sequences.
        # Default: False
        max_length = False,
        # Overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of
        # 1, i.e. even a single bp of overlapping sequence will be trimed off from the 3' end of any read.
        # Default: 1
        stringency = 1,
        # Maximum allowed error rate (no. of errors devided by the length of the matching region)
        # Default: 0.1
        e = 0.1,
        # Discard reads that became shorter than length INT because of either quality or adapter trimming.
        # A value of '0' effectively disables this behaviour.
        # For paired-end files, both reads of a read-pair need to be longer than <INT> bp to be printed out
        # to validated paired-end files (see option --paired).
        # If only one read became too short there is the possibility of keeping such unpaired single-end reads
        # (see --retain_unpaired).
        # Default: 20
        length = 20,
        # The total number of Ns (as integer) a read may contain before it will be removed altogether.
        # In a paired-end setting, either read exceeding this limit will result in the entire pair
        # being removed from the trimmed output files.
        # Default: False
        max_n = False,
        # Removes Ns from either side of the read. This option does currently not work in RRBS mode.
        # Default: False
        trim_n = False,
        # If specified no report file will be generated.
        # Default: False
        no_report_file = False,
        # If specified any output to STDOUT or STDERR will be wuppressed.
        # Default: False
        suppress_warn = False,
        # Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end reads).
        # This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at
        # the 5' end.
        # Default: False
        clip_R1 = False,
        # Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end reads only).
        # This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at
        # the 5' end. For paired-end BS-Seq, it is recommended to remove the first fewe bp because the end-repair
        # reaction may introduce a bias towards low methylation. Please refer to the M-bias plot section in the
        # Bismark User Guide for some examples.
        # Default: False
        clip_R2 = False,
        # Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end reads) AFTER
        # adapter/quality trimming has been performed. This may remove some unwanted bias from the 3' end
        # that is not directly related to adapter sequence or basecall quality.
        # Default: False
        three_prime_clip_R1 = False,
        # Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER adapter/quality trimming
        # has been performed. This may remove some unwanted bias from the 3' end that is not directly related
        # to adapter sequence or basecall quality.
        # Default: False
        three_prime_clip_R2 = False,
        # This enables the option '--nextseq-trim=3'CUTOFF' within Cutadapt, which will set a quality cutoff
        # (that is normally given with -q instead), but qualities of G bases are ignored.
        # This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any
        # signal are called as high-quality G bases. This is mutually exclusive with '-q INT'.
        nextseq = False,
        # Use PREFERRED_NAME as the basename for output files, instead of deriving the filenames from the
        # input files. Single-end data would be called PREFERRED_NAME_trimmed.fa(.gz), or
        # PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz) for paired-end data.
        # --basename only works when 1 file (single-end) or 2 files (paired-end) are specified, but nor for
        # longer lists.
        # Default: False
        basename = False,
        # 
        # RRBS-specific options (MspI digested material)
        # 
        # Specifies that the input file was an MspI digested RRBS sample (recognition site: CCGG).
        # Single-end or Read 1 sequences (paired-end) which were adapter-trimmed will have a further 2 bp 
        # removed from their 3' end. Sequences which were merely trimmed because of poor quality will not
        # be shortened further. Read 2 of paired-end libraries will in addition have the first 2 bp removed
        # from the 5' end (by setting '--clip_r2 2'). This is to avoid using artificial methylation calls from
        # the filled-in cytosine positions close to the 3' MspI site in sequenced fragments. This otion is
        # not recommended for users of the NuGEN ovation RRBS System 1-16 kit (see below).
        # Default: False
        rrbs = False,
        # Selecting this option for non-directional RRBS libraries will screen quality-trimmed sequences for
        # 'CAA' or 'CGA' at the start of the read and, if found, removes the first bwo basepairs. Like with
        # the option '--rrbs' this avoids using cytosine positions that were filled-in during the end-repair
        # step. '--non_directional' requires '--rrbs' to be specified as well. Note that this option does not
        # set '--clip_r2 2' in paired-end mode.
        # Default: False
        non_directional = False,
        # Keep the quality trimmed intermediate file, which means the temporary file is being deleted after
        # adapter trimming. Only has an effect for RRBS samples since other FastQ files are not trimmed for
        # poor qualities separately.
        # Default: False
        keep = False,
        #
        # Paired-end specific options
        #
        # Trims 1 bp off every read from its 3' end. This may be needed for FastQ files that are to be
        # aligned as paired-end data with Bowtie. This is because Bowtie regards alignments like this:
        #
        # R1 -------> or this:  ---------> R1
        # R2 <-------              <------ R2 
        # 
        # as invalid (whenever a start/end coordinate is contained within the other read).
        # NOTE: If you are planning to use Bowtie2, BWA etc. you don't need to specify this option.
        # Default: False
        trim1 = False,
        # If only one of the two paired-end reads became too short, the longer read will be written to
        # either '.unpaired_1.fq' or '.unpaired_2.fq' output files. The length cutoff for unpaired single-end
        # reads is governed by the parameters -r1/--length_1 and -r2/--length_2.
        # Default: False
        retain_unpaired = False,
        # Unpaired single-end read length cutoff needed for read 1 to be written to '.unpaired_1.fq' output
        # file. These reads may be mapped in single-end mode.
        # Default: 35
        length_1 = 35,
        # Unpaired single-end read length cutoff needed for read 2 to be written to '.unpaired_2.fq' output
        # file. These reads may be mapped in single-end mode.
        # Default: 35
        length_2 = 35,
    threads: 1
    log: 'logs/trim_galore/{sample}.log'
    benchmark: 'benchmarks/trim_galore/{sample}.log'
    wrapper: 'http://dohlee-bio.info:9193/trim-galore'

rule trim_galore_pe:
    input:
        '{sample}.read1.fastq.gz',
        '{sample}.read2.fastq.gz',
    output:
        '{sample}.read1.trimmed.fastq.gz',
        '{sample}.read2.trimmed.fastq.gz',
    params:
        extra = '',
        # Trim low-quality ends from reads in addition to adapter removal. For RRBS samples,
        # quality trimming will be performed first, and adapter trimming is carried in a second round.
        # Other files are quality and adapter trimmed in a single pass. The algorithm is the same
        # as the one used by BWA (Subtract INT from all qualities; compute partial sums from all
        # indices to the end of the sequence; cut sequence at the index at which the sum is minimal).
        # Default: 20
        quality = 20,
        # Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding)
        # for quality trimming.
        # Default: True
        phred33 = True,
        # Instructs Cutadapt to use ASCII+64 quality scores as Phred scores (Illumina 1.5 encoding) for
        # quality trimming.
        # Default: False
        phred64 = False,
        # Run FastQC in the default mode on the FastQ file once trimming is complete.
        # Default: False
        fastqc = False,
        # Passes extra arguments to FastQC. If more than one argument is to be passed to FastQC they must
        # be in the form "arg1 arg2 etc.". An example would be: --fastqc_args "--nogroup --outdir /home/".
        # Passing extra arguments will automatically invoke FastQC, so --fastqc does not have to be specified
        # separately.
        # Default: False
        fastq_args = False,
        # Adapter sequence to be trimmed. If no specified explicitly, Trim Galore will try to auto-detect
        # whether the Illumina universal, Nextera transposase or Illumina small RNA adapter sequence was used.
        # Also see '--illumina', '--nextera' and '--small_rna'. If no adapter can be detected within the
        # first 1 million sequences of the first file specified Trim Galore defaults to '--illumina'.
        # A single base may also be given as e.g. -a A{10}, to be expanded to -a AAAAAAAAAA.
        # Default: False
        adapter = False,
        # Optional adapter sequence to be trimmed off read 2 of paired-end files. This option requires
        # '--paired' to be specified as well. If the libraries to be trimmed are smallRNA then a2 will
        # be set to the Illumina small RNA 5' adapter automatically (GATCGTCGGACT). A single base may
        # also be given as e.g. -a2 A{10}, to be expanded to -a2 AAAAAAAAAA.
        # Default: False
        adapter2 = False,
        # Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter 'AGATCGGGAGC'
        # instead of the default auto-detection of adapter sequence.
        # Default: False
        illumina = False,
        # Adapter sequence to be trimmed is the first 12bp of the Nextera adapter 'CTGTCTCTTATA' instead of
        # default auto-detection of adapter sequence.
        # Default: False
        nextera = False,
        # Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter
        # 'TGGTTCTCGG' instead of the default auto-detection of adapter sequence. Selecting to trim
        # smallRNA adapters will also lower the --length value to 18bp. If the smallRNA libraries
        # are paired-end then a2 will be set to the Illumina small RNA 5' adapter automatically
        # (GATCGTCGGACT) unless -a 2 had been defined explicitly.
        # Default: False
        small_rna = False,
        # Discard reads that are longer than <INT> bp after trimming. This is only advised for smallRNA
        # sequencing to remove non-small RNA sequences.
        # Default: False
        max_length = False,
        # Overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of
        # 1, i.e. even a single bp of overlapping sequence will be trimed off from the 3' end of any read.
        # Default: 1
        stringency = 1,
        # Maximum allowed error rate (no. of errors devided by the length of the matching region)
        # Default: 0.1
        e = 0.1,
        # Discard reads that became shorter than length INT because of either quality or adapter trimming.
        # A value of '0' effectively disables this behaviour.
        # For paired-end files, both reads of a read-pair need to be longer than <INT> bp to be printed out
        # to validated paired-end files (see option --paired).
        # If only one read became too short there is the possibility of keeping such unpaired single-end reads
        # (see --retain_unpaired).
        # Default: 20
        length = 20,
        # The total number of Ns (as integer) a read may contain before it will be removed altogether.
        # In a paired-end setting, either read exceeding this limit will result in the entire pair
        # being removed from the trimmed output files.
        # Default: False
        max_n = False,
        # Removes Ns from either side of the read. This option does currently not work in RRBS mode.
        # Default: False
        trim_n = False,
        # If specified no report file will be generated.
        # Default: False
        no_report_file = False,
        # If specified any output to STDOUT or STDERR will be wuppressed.
        # Default: False
        suppress_warn = False,
        # Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end reads).
        # This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at
        # the 5' end.
        # Default: False
        clip_R1 = False,
        # Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end reads only).
        # This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at
        # the 5' end. For paired-end BS-Seq, it is recommended to remove the first fewe bp because the end-repair
        # reaction may introduce a bias towards low methylation. Please refer to the M-bias plot section in the
        # Bismark User Guide for some examples.
        # Default: False
        clip_R2 = False,
        # Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end reads) AFTER
        # adapter/quality trimming has been performed. This may remove some unwanted bias from the 3' end
        # that is not directly related to adapter sequence or basecall quality.
        # Default: False
        three_prime_clip_R1 = False,
        # Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER adapter/quality trimming
        # has been performed. This may remove some unwanted bias from the 3' end that is not directly related
        # to adapter sequence or basecall quality.
        # Default: False
        three_prime_clip_R2 = False,
        # This enables the option '--nextseq-trim=3'CUTOFF' within Cutadapt, which will set a quality cutoff
        # (that is normally given with -q instead), but qualities of G bases are ignored.
        # This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any
        # signal are called as high-quality G bases. This is mutually exclusive with '-q INT'.
        nextseq = False,
        # Use PREFERRED_NAME as the basename for output files, instead of deriving the filenames from the
        # input files. Single-end data would be called PREFERRED_NAME_trimmed.fa(.gz), or
        # PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz) for paired-end data.
        # --basename only works when 1 file (single-end) or 2 files (paired-end) are specified, but nor for
        # longer lists.
        # Default: False
        basename = False,
        # 
        # RRBS-specific options (MspI digested material)
        # 
        # Specifies that the input file was an MspI digested RRBS sample (recognition site: CCGG).
        # Single-end or Read 1 sequences (paired-end) which were adapter-trimmed will have a further 2 bp 
        # removed from their 3' end. Sequences which were merely trimmed because of poor quality will not
        # be shortened further. Read 2 of paired-end libraries will in addition have the first 2 bp removed
        # from the 5' end (by setting '--clip_r2 2'). This is to avoid using artificial methylation calls from
        # the filled-in cytosine positions close to the 3' MspI site in sequenced fragments. This otion is
        # not recommended for users of the NuGEN ovation RRBS System 1-16 kit (see below).
        # Default: False
        rrbs = False,
        # Selecting this option for non-directional RRBS libraries will screen quality-trimmed sequences for
        # 'CAA' or 'CGA' at the start of the read and, if found, removes the first bwo basepairs. Like with
        # the option '--rrbs' this avoids using cytosine positions that were filled-in during the end-repair
        # step. '--non_directional' requires '--rrbs' to be specified as well. Note that this option does not
        # set '--clip_r2 2' in paired-end mode.
        # Default: False
        non_directional = False,
        # Keep the quality trimmed intermediate file, which means the temporary file is being deleted after
        # adapter trimming. Only has an effect for RRBS samples since other FastQ files are not trimmed for
        # poor qualities separately.
        # Default: False
        keep = False,
        #
        # Paired-end specific options
        #
        # This option performs length trimmng of quality/adapter/RRBS trimmed reads for paired-end files. To
        # pass the validation test, both sequences of a sequence pair are required to have a certain minimum
        # length which is governed by the option --length (see above). If only one read passes this length
        # threshold the other read cen be rescued (see option --retain_unpaired). Using this option lets 
        # you discard too short read pairs without disturbint the sequence-by-sequence order of FastQ files
        # which is required by many aligners.
        # Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g. file1_1.fq
        # file1_2.fq, SRR2_1.fa.gz SRR2_2.fq.gz, ...
        # Default: False
        paired = False,
        # Trims 1 bp off every read from its 3' end. This may be needed for FastQ files that are to be
        # aligned as paired-end data with Bowtie. This is because Bowtie regards alignments like this:
        #
        # R1 -------> or this:  ---------> R1
        # R2 <-------              <------ R2 
        # 
        # as invalid (whenever a start/end coordinate is contained within the other read).
        # NOTE: If you are planning to use Bowtie2, BWA etc. you don't need to specify this option.
        # Default: False
        trim1 = False,
        # If only one of the two paired-end reads became too short, the longer read will be written to
        # either '.unpaired_1.fq' or '.unpaired_2.fq' output files. The length cutoff for unpaired single-end
        # reads is governed by the parameters -r1/--length_1 and -r2/--length_2.
        # Default: False
        retain_unpaired = False,
        # Unpaired single-end read length cutoff needed for read 1 to be written to '.unpaired_1.fq' output
        # file. These reads may be mapped in single-end mode.
        # Default: 35
        length_1 = 35,
        # Unpaired single-end read length cutoff needed for read 2 to be written to '.unpaired_2.fq' output
        # file. These reads may be mapped in single-end mode.
        # Default: 35
        length_2 = 35,
    threads: 1
    log: 'logs/trim_galore/{sample}.log'
    benchmark: 'benchmarks/trim_galore/{sample}.log'
    wrapper: 'http://dohlee-bio.info:9193/trim-galore'
