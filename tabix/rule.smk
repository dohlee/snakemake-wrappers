rule tabix_index_vcf:
    # NOTE: When executed, this rule will overwrite existing index without asking.
    input:
        '{file}.vcf.gz'
    output:
        '{file}.vcf.gz.tbi'
    params:
        extra = ''
    threads: 1
    logs: 'logs/tabix/{file}.log'

rule tabix_index_bed:
    # NOTE: When executed, this rule will overwrite existing index without asking.
    input:
        '{file}.bed.gz'
    output:
        '{file}.bed.gz.tbi'
    params:
        extra = ''
    threads: 1
    logs: 'logs/tabix/{file}.log'

rule tabix_index_gff:
    # NOTE: When executed, this rule will overwrite existing index without asking.
    input:
        '{file}.gff.gz'
    output:
        '{file}.gff.gz.tbi'
    params:
        extra = ''
    threads: 1
    logs: 'logs/tabix/{file}.log'
