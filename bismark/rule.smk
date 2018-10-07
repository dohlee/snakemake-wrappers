def get_bismark_inputs(wildcards):
    """Define your function to tell whether the sample is
    single-ended or paired-ended.
    """
    raise NotImplementedError('Function get_bismark_inputs is not implemented.')

    # IMPLEMENT FUNCTION LIKE:
    if is_paired(s):
        return dict(
            fastq='data/{sample}.trimmed.fastq.gz',
            reference_dir=directory('reference/'),
            bisulfite_genome_dir=directory('reference/Bisulfite_Genome'),
        )
    else:
        return dict(
            fastq=[
                'data/{sample}.read1.trimmed.fastq',
                'data/{sample}.read2.trimmed.fastq',
            ],
            reference_dir=directory('reference/'),
            bisulfite_genome_dir=directory('refernce/Bisulfite_Genome'),
        )


rule bismark:
    input: unpack(get_bismark_inputs)
    output:
        'result/{sample}/{sample}.bismark.bam',
        'result/{sample}/{sample}.bismark_report.txt'
    threads: 6
    wrapper:
        'http://dohlee-bio.info:9193/bismark'
