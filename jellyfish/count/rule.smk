rule jellyfish_count:
    input:
        # Required input.
        fastq = ['data/{sample}.read1.fastq.gz', 'data/{sample}.read2.fastq.gz']
    output:
        # Required output.
        'result/{sample}.jf',
    params:
        extra = '',
        # (Required) Length of mer.
        mer_len = 31,
        # (Required) Initial hash size.
        size = 800000000,
        # SAM/BAM/CRAM formatted input file.
        sam = False,
        # Number files open simultaneously.
        Files = False,
        # File of commands generating fast[aq].
        generator = False,
        # Number of generators run simultaneously.
        Generators = False,
        # Shell used to run generator commands.
        shell_ = False,
        # Length bits of counting field.
        # Default: 7
        counter_len = 7,
        # Length in bytes of counter field in output.
        # Default: 4
        out_counter_len = 4,
        # Count both strand, canonical representation.
        canonical = False,
        # Bloom counter to filter out singleton mers.
        bc = False,
        # Use bloom filter to count high-frequency mers.
        bf_size = False,
        # False positive rate of bloom filter.
        # Default: 0.01
        bf_fp = False,
        # Count only k-mers in this files.
        if_ = False,
        # Any base with quality below this character is changed to N.
        min_qual_char = False,
        # ASCII for quality values.
        # Default: 64
        quality_start = 64,
        # Minimum quality. A base with lesser quality becomes an N.
        min_quality = False,
        # Maximum number of reprobes.
        # Default: 126
        reprobes = 126,
        # Dump in text format.
        text = False,
        # Disk operation. Do not do size doubling.
        disk = False,
        # Don't output k-mer with count < lower-count.
        lower_count = False,
        # Don't output k-mer with count > upper-count.
        upper_count = False,
        # Print timing information.
        timing = False,
    threads: 8
    log: 'logs/jellyfish_count/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/jellyfish/count'
