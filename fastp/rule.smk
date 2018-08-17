rule fastp_se:
    input:
        '{sample}.fastq.gz'
    output:
        '{sample}.fastp_filtered.fastq.gz'
    params:
        extra = '' \
        # Reporting options
        # The json format report file name.
        # '-j STR.json ' \
        # The html format report file name.
        # '-h STR.html ' \
        # Set compression level.
        # '-Z INT ' \
        # Disable adapter trimming.
        # '-A ' \
        # The adapter for r1.
        # '--adapter_sequence ' \
        # The adapter for r2.
        # '--adapter_sequence_r2 ' \
        # Disable quality filtering.
        # '-Q '
    wrapper:
        'http://dohlee-bio.info:9193/fastp'

rule fastp_pe:
    input:
        '{sample}.read1.fastq.gz',
        '{sample}.read2.fastq.gz'
    output:
        '{sample}.read1.fastp_filtered.fastq.gz',
        '{sample}.read2.fastp_filtered.fastq.gz'
    params:
        extra = '' \
        # Reporting options
        # The json format report file name.
        # '-j STR.json ' \
        # The html format report file name.
        # '-h STR.html ' \
        # Set compression level.
        # '-Z INT ' \
        # Disable adapter trimming.
        # '-A ' \
        # The adapter for r1.
        # '--adapter_sequence ' \
        # The adapter for r2.
        # '--adapter_sequence_r2 ' \
        # Disable quality filtering.
        # '-Q '
    wrapper:
        'http://dohlee-bio.info:9193/fastp'
