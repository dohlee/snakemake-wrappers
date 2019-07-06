rule wxs_normal:
    output: ['normal.read1.fastq.gz', 'normal.read2.fastq.gz']
    wrapper: 'http://dohlee-bio.info:9193/test/wxs/normal'
