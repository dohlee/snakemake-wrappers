rule reference:
    output:
        [
            'HBR1.read1.fastq.gz', 'HBR1.read2.fastq.gz',
            'HBR2.read1.fastq.gz', 'HBR2.read2.fastq.gz',
            'HBR3.read1.fastq.gz', 'HBR3.read2.fastq.gz',
            'UHR1.read1.fastq.gz', 'UHR1.read2.fastq.gz',
            'UHR2.read1.fastq.gz', 'UHR2.read2.fastq.gz',
            'UHR3.read1.fastq.gz', 'UHR3.read2.fastq.gz',
        ]
    wrapper:
        'http://dohlee-bio.info/test/rna-seq/pe/deg'
