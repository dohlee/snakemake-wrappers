rule bwa_mem:
    input:
        # Required input. Input read file.
        reads = ['{sample}.read1.fastq.gz', '{sample}.read2.fastq.gz'],
        # You may use any of {genome}.amb, {genome}.ann, {genome}.bwt,
        # {genome}.pac, {genome}.sa just to obligate snakemake to run `bwa index` first.
        reference = 'reference/hg38.bwt'
    output:
        # BAM output or SAM output is both allowed.
        # Note that BAM output will be automatically detected by its file extension,
        # and SAM output (which is bwa mem default) will be piped through `samtools view`
        # to convert SAM to BAM.
        'result/{sample}/{sample}.bam'
    params:
        extra = '',
        # Minimum seed length.
        # Default: 19
        k = 19,
        # Band width for banded alignment.
        # Default: 100
        w = 100,
        # Off-diagonal X-dropoff.
        # Default: 100
        d = 100,
        # Look for internal seeds inside a seed longer than {-k} * FLOAT
        # Default: 1.5
        r = 1.5,
        # Seed occurrence for the 3rd round seeding.
        # Default: 20
        y = 20,
        # Skip seeds with more than INT occurrences.
        # Default: 500
        c = 500,
        # Drop chains shorter than FLOAT fraction of the logest overlapping chain.
        # Default: 0.5
        D = 0.50,
        # Discard a chain if seeded bases shorter than INT.
        # Default: 0
        W = 0,
        # Perform at most INT rounds of mate rescues for each read.
        # Default: 50
        m = 50,
        # Skip mate rescue.
        # Default: False
        S = False,
        # Skip pairing; mate rescue performed unless -S also in use
        # Default: False
        P = False,
        # Score for a sequence match, which scales options -TdBOELU unless overridden.
        # Default: 1
        A = 1,
        # Penalty for a mismatch.
        # Default: 4
        B = 4,
        # Gap open penalties for deletions and insertions.
        # Default: 6,6
        O = '6,6',
        # Gap extension penalty; a gap of size k cost '{-O} + {-E}*k'.
        # Default: 1,1
        E = '1,1',
        # Penalty for 5'- and 3'-end clipping.
        # Default: 5,5
        L = '5,5',
        # Penalty for an unpaired read pair.
        # Default: 17
        U = 17,
        # Read type. Setting -x changes multiple parameters unless overridden [null]
        # pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 (PacBio reads to ref)
        # ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0 (Oxford Nanopore 2D-reads to ref)
        # intractg: -B9 -O16 -L5 (intra-species contigs to ref)
        # Default: False
        x = False,
        # Smart pairing (ignoring in2.fq)
        # Default: False
        p = False,
        # Read group header line such as '@RG\tID:foo\tSM:bar'
        # Default: False
        # NOTE: You should check the platform information of the read data!
        R = r"'@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        # Insert STR to header if it starts with @; or insert lines in FILE.
        # Default: False
        H = False,
        # Treat ALT contigs as part of the primary assembly. (i.e. ignore <idxbase>.alt file)
        # Default: False
        j = False,
        # For split alignment, take the alignment with the smallest coordinate as primary.
        # Default: False
        _5 = False,
        # Dont't modify mapQ of supplementary alignments.
        # Default: False
        q = False,
        # Process INT input bases in each batch regardless of nThreads (for reproducibility).
        # Default: False.
        K = False,
        # Verbosity level: 1=error, 2=warning, 3=message, 4+=debugging
        # Default: 3
        v = 3,
        # Minimum score to output.
        # Default: 30
        T = 30,
        # If there are <INT hits with score > 80% of the max score, output all in XA.
        # Default: 5,200
        h = '5,200',
        # Output all alignments for SE or unpaired PE.
        # Default: False
        a = False,
        # Append FASTA/FASTQ comment to SAM output.
        # Default: False
        C = False,
        # Output the reference FASTA header in the XR tag.
        # Default: False
        V = False,
        # Use soft clipping for supplementary alignments.
        # Default: False
        Y = False,
        # Mark shorter split hits as secondary.
        # NOTE: You may need this if you use GATK downstream.
        # Default: False
        M = False,
        # Specify the mean, standard deviation (10% of the mean if absent), max
        # (4 sigma from the mean if absent) and min of the insert size distribution.
        # FR orientation only.
        # Default: False (inferred)
        I = False,
    threads: 8
    log: 'logs/bwa_mem/{sample}.log'
    benchmark: 'benchmarks/bwa_mem/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/bwa/mem'
