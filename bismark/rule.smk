rule bismark_single:
    input:
        fastq = 'data/{sample}.fastq.gz',
        reference_dir = directory('reference/'),
        bisulfite_genome_dir = directory('reference/Bisulfite_Genome')
    output:
        'result/{sample}/{sample}_bismark_bt2.bam',
        'result/{sample}/{sample}_bismark_bt2_SE_report.txt'
    threads: 6
    wrapper:
        'http://dohlee-bio.info:9193/bismark'

rule bismark_paired:
    input:
        fastq = ['data/{sample}.read1.fastq', 'data/{sample}.read2.fastq'],
        reference_dir = directory('reference/'),
        bisulfite_genome_dir = directory('reference/Bisulfite_Genome')
    output:
        'result/{sample}/{sample}_bismark_bt2_pe.bam',
        'result/{sample}/{sample}_bismark_bt2_PE_report.txt'
    wrapper:
        'http://dohlee-bio.info:9193/bismark'
