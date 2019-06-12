rule bismark_genome_preparation:
    input:
        directory('reference/')
    output:
        directory('reference/Bisulfite_Genome')
    params:
        extra = '',
        # Print verbose output for more details or debugging.
        # Default: False
        verbose = False,
        # The full path to the Bowtie2 or HISAT2 installation folder on your system (depending on which
        # aligner/indexer you intend to use; please note that this is the folder and not any executable).
        # Unless this path is specified, it is assumed that the aligner in question (Bowtie2/HISAT2) is
        # in the PATH.
        # Default: False
        path_to_aligner = False,
        # This will create bisulfite indexes for use with Bowtie2. Recommended for most bisulfite sequencing
        # applications.
        # Default: True
        bowtie2 = True,
        # This will create bisulfite indexes for use with HISAT2. At the time of writing, this is
        # still unchartered territory, and only reommended for specialist applications such as RNA-methylation
        # analyses or SLAM-seq type applications (see also: --slam).
        # Default: False
        hisat2 = False,
        # Instruct the Bismark Indexer to write the converted genomes into single-entry FastA files instead
        # of making one multi-FastA file (MFA) per chromosome. This might be useful if individual bisulfite
        # converted chromosomes are needed (e.g. for debugging), however it can cause a problem with indexing
        # if the number of chromosomes is vast (this is likely to be in the range of several thousand files;
        # the operating system can only handle lists up to a ceratin length, and some newly assembled genomes
        # may contain 20000-500000 contigs of scaffold files which do exceed this list length limit).
        # Default: False
        single_fasta = False,
        # Calculate and extract the genomic sequence composition for mono and di-nucleotides and write the
        # genomic composition table 'genomic_nucleotide_frequencies.txt' to the genome folder. This may be
        # useful later on when using bam2nuc or the Bismark option --nucleotide_coverage
        # Default: False
        genomic_composition = False,
        # Instead of performing an in-silico bisulfite conversion, this mode transforms T to C (forward 
        # strand), or A to G (reverse strand). The folder structure and rest of the indexing prcess is
        # currently exactly the smae as for bisulfite sequences, but this might change at some point.
        # This means that a genome prepared in --slam mode is currently indistinguishable from a true
        # Bisulfite Genome, so please make sure you name the genome folder appropriately to avoid confusion.
        # Default: False
        slam = False,
        # Force generated index to be 'arge', even if reference has fewer than 4 billion nucleotides.
        # At the time of writing this is required for parallel processing of VERY LARGE genomes (e.g.
        # the axolotl).
        # Default: False
        large_index = False,
    threads: 1
    log: 'logs/bismark_genome_preparation/_.log'
    benchmark: 'benchmarks/bismark_genome_prepration/_.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/bismark/genome-preparation'
