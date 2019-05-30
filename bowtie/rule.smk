rule bowtie:
    input:
        index_dir = 'DIRECTORY_FOR_YOUR_BOWTIE_INDEX',
        # For single reads,
        #reads = [
        #    '{sample}.fastq.gz',
        #],
        # For paired-end reads,
        reads = [
            '{sample}.read1.fastq.gz',
            '{sample}.read2.fastq.gz',
        ]
    output:
        # It automatically sorts the output bam file if its file name ends with '.sorted.bam',
        # e.g.
        # '{sample}.sorted.bam'
        bam = '{sample}.bam'
    params:
        # Additional parameters go here.
        extra = '',
        # Seed for random number generator.
        seed = 0,
        # Skip the first <int> reads/pairs in the input.
        skip = False,
        # Stop after first <int> reads/pairs, excluding skipped reads.
        qupto = False,
        # Trim <int> bases from 5' (left) end of reads
        trim5 = False,
        # Trim <int> bases from 3' (right) end of reads
        trim3 = False,
        # Input quals are Phred+33 (default)
        phred33_quals = True,
        # Input quals are Phred+64 (same as --solexa1.3-quals)
        phred64_quals = False,
        # Input quals are from GA Pipeline ver. < 1.3
        solexa_quals = False,
        # Input quals are from GA Pipeline ver. >= 1.3
        solexa13_quals = False,
        # Qualities are given as space-separated integers (not ASCII)
        integer_quals = False,
        # Force usage of a 'arge' index, even if a small one is present.
        large_index = False,
        # Report end-to-end hits with <=v mismatches; ignore qualities.
        v = False,
        # Max mismatches in seed (can be 0-3)
        # Deafult: 2
        seedmms = 2,
        # Max sum of mismatch quals across alignment for -n
        # Default: 70
        maqerr = 70,
        # Seed length for -n
        # Default: 28
        seedlen = 28,
        # Disable Maq-like quality rounding for -n (nearest 10 <= 30)
        nomaqround = False,
        # Minimum insert size for paired-end alignment
        # Default: 0
        minins = 0,
        # Maximum insert size for paired-end alignment
        # Default: 250
        maxins = 250,
        # -1, -2 mates align fw/rev, rev/fw, fw/fw.
        # Default: fr
        fr = True,
        rf = False,
        ff = False,
        # Do not align to forward/reverse-complement reference strand.
        nofw = False,
        norc = False,
        # Max # backtracks for -n 2/3
        # Default: 125, 800 for --best
        maxbts = False,
        # Max # attempts to find mate for anchor hit.
        # Default: 100
        pairtries = 100,
        # Try hard to find valid alignments, at the expense of speed.
        tryhard = False,
        # Max megabytes of RAM for best-first search frames.
        # Default: 64
        chunkmbs = 64,
        # # of reads to read from input file at once.
        # Default: 16
        reads_per_batch = 16,
        # Report up to <int> good alignments per read.
        # Default: 1
        k = 1,
        # Report all alignment sper read (much slower than low -k)
        all = False,
        # Suppress all alignments if > <int> exist.
        # Default: no limit
        m = False,
        # Like -m, but reports 1 random hit (MAPQ=0); requires --best
        M = False,
        # Hist guaranteed best stratu; ties broken by quality.
        best = False,
        # Hits in sub-optimal strata arent' reported (requires --best)
        strata = False,
        # Print wall-clock time taken by search phases.
        time = False,
        # Leftmost ref offset = <int> in bowtie output.
        # Default: 0
        offbase = 0,
        # Print nothing but the alignments.
        quiet = False,
        # Refer to ref. seqs by 0-based index rather than name.
        refidx = False,
        # Write aligned reads/pairs to file(s) <fname>
        al = False,
        # Write unaligned reads/pairs to file(s) <fname>
        un = False,
        # Suppress SAM records for unaligned reads.
        no_unal = False,
        # Write reads/pairs over -m limit to file(s) <fname>
        max = False,
        # Suppresses given columns (comma-delim'ed) in default output.
        suppress = False,
        # Write entire ref name.
        # Default: Only up to 1st space
        fullref = False,
        # Phred penalty for SNP when decoding colorspace.
        # Default: 30
        snpphred = 30,
        # Approximate fraction of SNP bases (e.g. 0.001); sets --snpphred
        snpfrac = False,
        # Print aligned colorspace seqs as colors, not decoded bases.
        col_cseq = False,
        # Print original colorspace quals, not decoded quals.
        col_cqual = False,
        # Keep nucleotides at extreme ends of decoded alignment.
        col_keepends = False,
        # Write hits in SAM format.
        # Default: True (for this wrapper)
        sam = True,
        # Default mapping quality (MAPQ) to print for SAM alignments.
        mapq = False,
        # Suppress header lines (starting with @) for SAM output.
        sam_nohead = False,
        # Suppress @SQ header lines for SAM output.
        sam_nosq = False,
        # ADD <text> (usually "lab=value") to @RG line of SAM header.
        sam_rg = False,
        # Override offrate of index; must be >= index's offrate.
        offrate = False,
        # Use memory-mapped I/O for index; many 'bowtie's can share.
        mm = False,
        # Use shared mem for index; many 'bowtie's can share.
        shmem = False,
        # Discard mapped reads having mapping quality (MAPQ) below this value.
        # NOTE: This will be done after the alignment, with `samtools view -bS -q <int>` command.
        mapq_cutoff = 10,
    threads: 4
    benchmark:
        repeat('benchmarks/bowtie/{sample}.tsv', 1)
    log: 'logs/bowtie/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bowtie'

