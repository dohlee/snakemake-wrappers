rule wgsim:
    input:
        # Required input.
        'reference/Homo_sapiens.GRCh38.genome.fasta'
    output:
        # Required output.
        # You can specify single-read output, too.
        'result/simulated.read1.fastq',
        'result/simulated.read2.fastq'
    params:
        # Optional parameters. It omitted, default value will be used.
        random_seed = 1,  # Seed for random generator. [default: 0, which means to use the current time]
        N = 10000,  # Number of read pairs. [default: 1000000]
        error_rate = 0.01,  # Base error rate. [default: 0.020]
        mutation_rate = 0.001,  # Rate of mutations. [default: 0.001]
        read1_length = 70,  # Length of the first read. [default: 70]
        indel_fraction = 0.15,  # Fraction of indels. [default: 0.15]
        indel_extension_probability = 0.30,  # Probability an indel is extended. [default: 0.30]
        # Discard if the fraction of ambiguous bases higher than that. [default: 0.05]
        discard_ambiguous_bases_fraction = 0.05,
        haplotype_mode = False,
        extra = '',
        #
        # Paired-end only options.
        #
        read2_length = 70,  # Length of the second read. [default: 70]
        read_distance = 500,  # Outer distance between the two ends. [default: 500]
        standard_deviation = 50,  # Standard deviation. [default: 50]
    threads: 1
    log: 'logs/wgsim/wgsim.log'
    wrapper:
        'http://dohlee-bio.info:9193/wgsim'
