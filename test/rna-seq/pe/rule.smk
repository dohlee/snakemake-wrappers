rule reference:
    output:
        ['test_paired.read1.fastq.gz', 'test_paired.read2.fastq.gz']
    wrapper:
        'http://dohlee-bio.info/test/rna-seq/pe'
