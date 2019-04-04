from os.path import join

rule bowtie_build:
    input:
        # Required input.
        # Reference genome fasta.
        # e.g.
        # reference = config['reference']['genome']['fasta']
        reference = 'REFERENCE_FASTA'
    output:
        # e.g.
        # reference_indices = config['reference']['genome']['bowtie_index_dir']
        index_dir = directory('DIRECTORY_FOR_YOUR_BOWTIE_INDEX')
    params:
        # Additional parameters go here.
        extra = ''
    threads: 4
    benchmark:
        repeat("benchmarks/bowtie_build/{genome}.tsv", 1)
    log: 'logs/bowtie_build/{genome}.log'
    wrapper:
    'http://dohlee-bio.info:9193/bowtie/build'