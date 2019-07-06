rule wxs_tumor:
    output: ['tumor.read1.fastq.gz', 'tumor.read2.fastq.gz']
    wrapper: 'http://dohlee-bio.info:9193/test/wxs/tumor'
