rule rrbs_pe:
    output: ['rrbs_pe.read1.fastq.gz', 'rrbs_pe.read2.fastq.gz']
    wrapper: 'http://dohlee-bio.info:9193/test/rrbs/pe'
