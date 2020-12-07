from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'simulated.fastq.gz', 'simulated.read1.fastq.gz', 'simulated.read2.fastq.gz'

rule gzip_fastq:
    input: '{file}.fastq'
    output: '{file}.fastq.gz'
    shell: 'gzip {input}'

rule gunzip_fastq:
    input: '{file}.fastq.gz'
    output: '{file}.fastq'
    shell: 'gunzip {input}'

rule gunzip_fasta:
    input: '{file}.fasta.gz'
    output: '{file}.fasta'
    shell: 'gunzip {input}'

rule download_reference_genome:
    input: HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/reference/Homo_sapiens.GRCh38.genome.chr20.fasta.gz')
    output: 'Homo_sapiens.GRCh38.genome.chr20.fasta.gz'
    shell: 'mv {input} {output}'

rule wgsim_se:
    input:
        'Homo_sapiens.GRCh38.genome.chr20.fasta'
    output:
        'simulated.fastq'
    params:
        # Optional parameters. It omitted, default value will be used.
        random_seed = 1,  # Seed for random generator. [default: 0, which means to use the current time]
        N = 10000,  # Number of read pairs. [default: 1000000]
        error_rate = 0.01,  # Base error rate. [default: 0.020]
        mutation_rate = 0.001,  # Rate of mutations. [default: 0.001]
        standard_deviation = 50,  # Standard deviation. [default: 50]
        read1_length = 70,  # Length of the first read. [default: 70]
        indel_fraction = 0.15,  # Fraction of indels. [default: 0.15]
        indel_extension_probability = 0.30,  # Probability an indel is extended. [default: 0.30]
        # Discard if the fraction of ambiguous bases higher than that. [default: 0.05]
        discard_ambiguous_bases_fraction = 0.05,
        haplotype_mode = False,
        extra = '',
    threads: 1
    log:
        'logs/wgsim/wgsim_se.log'
    wrapper:
        'wgsim'

rule wgsim_pe:
    input:
        'Homo_sapiens.GRCh38.genome.chr20.fasta'
    output:
        'simulated.read1.fastq',
        'simulated.read2.fastq'
    params:
        # Optional parameters. It omitted, default value will be used.
        random_seed = 1,  # Seed for random generator. [default: 0, which means to use the current time]
        N = 10000,  # Number of read pairs. [default: 1000000]
        error_rate = 0.01,  # Base error rate. [default: 0.020]
        mutation_rate = 0.001,  # Rate of mutations. [default: 0.001]
        standard_deviation = 50,  # Standard deviation. [default: 50]
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
    threads: 1
    log:
        'logs/wgsim/wgsim_pe.log'
    wrapper:
        'wgsim'
