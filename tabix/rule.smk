rule tabix_index:
    # NOTE: When executed, this rule will overwrite existing index without asking.
    # input should end with either of '.vcf.gz', '.bed.gz', or '.gff.gz'.
    input:
        '{file}.vcf.gz'
        # '{file}.bed.gz'
        # '{file}.gff.gz'
    output:
        '{file}.vcf.gz.tbi'
        # '{file}.bed.gz.tbi'
        # '{file}.gff.gz.tbi'
    params:
        extra = ''
    threads: 1
    logs: 'logs/tabix/{file}.log'
    wrapper: 'http://dohlee-bio.info:9193/tabix'