from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: 'result/test/abundance.tsv'

rule fetch_reference_transcriptome:
    input: HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/reference/Homo_sapiens.GRCh38.transcriptome.chr20.fasta.gz')
    output: 'Homo_sapiens.GRCh38.transcriptome.chr20.fasta'
    shell: 'gunzip -c {input} > {output}'

rule fetch_rna_seq_se_fastq:
    input: HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/rna-seq/se/test.fastq.gz')
    output: 'test.fastq.gz'
    shell: 'mv {input} {output}'

rule kallisto_index:
    input:
        # Required input.
        'Homo_sapiens.GRCh38.transcriptome.chr20.fasta'
    output:
        # Required output.
        'Homo_sapiens.GRCh38.transcriptome.chr20.kallisto_idx'
    params:
        # Optional parameters. It omitted, default value will be used.
        extra = '',
    threads: 1
    log: 'logs/kallisto_index/Homo_sapiens.GRCh38.genome.log'
    wrapper:
        'kallisto/index'

rule kallisto_quant:
    input:
        # Required input.
        index = 'Homo_sapiens.GRCh38.transcriptome.chr20.kallisto_idx',
        fastq = ['{sample}.fastq.gz']
    output:
        # Required output.
        'result/{sample}/abundance.tsv',
        'result/{sample}/abundance.h5',
        'result/{sample}/run_info.json'
    params:
        fragment_length = 75.0,
        standard_deviation = 10.0,
        # Optional parameters. It omitted, default value will be used.
        extra = '',
    threads: 1
    log: 'logs/kallisto_quant/{sample}.log'
    wrapper:
        'kallisto/quant'
