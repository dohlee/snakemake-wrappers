rule rsem_calculate_expression:
    input:
        # For paired-end reads, use:
        reads = ['{sample}.read1.fastq.gz', '{sample}.read2.fastq.gz'],
        # For single-end reads, use:
        # reas = ['{sample}.fastq.gz'],
        # NOTE: REFERENCE_PREFIX.transcipts.fa should be generated via
        # rsem-prepare-reference
        reference = 'REFERENCE_PREFIX.transcripts.fa',
    output:
        genes = '{sample}.genes.results',
        isoforms = '{sample}.isoforms.results'
    params:
        extra = lambda wildcards: 'PAIRED' if SAMPLE2LIB[wildcards.sample] == 'PAIRED' else '',
        # Input reads do not contain quality scores.
        # If this option is on, fasta input is assumed.
        no_qualities = False,
        # This option defines the strandedness of the RNA-Seq reads.
        # It recognizes three values: 'none', 'forward', and 'reverse'.
        # 'none': Non-strand-specific protocols.
        # 'forward': All (upstream) reads are derived from the forward strand.
        # 'reverse': All (upstream) reads are derived from the reverse strand.
        # If forward/reverse is set, the '--norc'/'--nofw' Bowtie/Bowtie2 option
        # will also be enabled to avoid aligning reads to the opposite strand.
        # For Illumina truSeq Stranded protocols, please use 'reverse'.
        # Default: 'none'
        strandedness = 'none',
        # Input file contains alignments in SAM/BAM/CRAM format.
        # The exact file format will be determined automatically.
        # Default: False
        alignments = False,
        # Use Bowtie2 instead of Bowtie to align reads.
        # Default: False
        bowtie2 = False,
        # Use STAR to align reads. Alignment parameters are from ENCODE3's
        # STAR-RSEM pipeline.
        # Default: False
        star = False,
        # If gene_name/transcript_name is available, append it to the end of
        # gene_id/transcript_id (separated by '_') in files
        # 'sample_name.isoforms.results' and 'sample_name.genes.results'
        # Default: False
        append_names = False,
        # Set the seed for the random number generators used in calculating
        # posterior mean estimates and credibility intervals.
        # The seed must be a non-negative 32 bit integer.
        # Default: 0
        seed = 0,
        # By default, RSEM uses Dirichlet(1) as the prior to calculate
        # posterior mean estimates and credibility intervals. However, much less
        # genes are expressed in single cell RNA-Seq data. Thus, if you want
        # to compute posterior mean estimates and/or credibility intervals
        # and you have single-cell RNA-Seq data, you are recommended to turn
        # this option. Then RSEM will use Dirichlet(0.1) as the prior which
        # encourage the sparsity of the expression levels.
        # Default: False
        single_cell_prior = False,
        # Run RSEM's collapsed Gibbs sampler to calculate posterior mean estimates.
        # Default: False
        calc_pme = False,
        # Calculate 95% credibility intervals and posterior mean esimates.
        # The credibility level can be changed by setting '--ci-credibility-level'.
        # Default: False
        calc_ci = False,
        # Supress the output of logging information.
        # Default: False
        quiet = False,
        # Sort BAM file aligned under transcript coordinate by read name.
        # Setting this option on will produce deterministic maximum likelihood
        # estimations from independent runs. Note that osorting will take long
        # time and lots of memory.
        # Default: False
        sort_bam_by_read_name = False,
        # Do not output any BAM file.
        # Default: False
        no_bam_output = False,
        # When RSEM generates a BAM file, instead of outputting all alignments a
        # read has with their posterior probabilities, one alignment is sampled
        # according to the posterior probabilities. The sampling procedure
        # includes the alignment to the "noise" transcript, which does not
        # appear in the BAM file. Only the sampled alignment has a weight of 1.
        # All other alignments have weight 0. If the "noise" transcript is
        # sampled, all alignments appeared in the BAM file should have weight 0.
        # Default: False
        sampling_for_bam = False
        # Generate a BAM file, 'sample_name.genome.bam', with alignments ampped
        # to genoic coordinates and annotated with their posterior probabilities.
        # In addition, RSEM will call samtools (included in RSEM package) to sort
        # and index the bam file.
        # 'sample_name.genome.sorted.bam' and 'sample_name.genome.sorted.bam.bai'
        # will be generated.
        # Default: False
        output_genome_bam = False,
        # Sort RSEM generated transcript and genome BAM files by coordinates and
        # build associated indices.
        # Default: False
        sort_bam_by_coordinate = False,
        # Set the maximum memory per thread that can be used by 'samtools sort',
        # Maximum memory accepts sufficies 'K/M/G'. RSEM will pass it to the '-m'
        # option of 'samtools sort'. Note that the default used here is different
        # from the default used by samtools.
        # Default '1G'.
        sort_bam_memory_per_thread = '1G',
    threads: 8
    log: 'logs/rsem_calculate_expression/{sample}.log'
    benchmark: 'benchmarks/rsem_calculate_expression/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/rsem/calculate-expression'
