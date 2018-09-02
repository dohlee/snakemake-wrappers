rule hisat2_build:
    # WARNING: If you use SNP, haplotype, splice site, or exon information,
    # hisat2-build requires 200GB of RAM for indexing.
    # If that much of RAM is not available, you can download
    # pre-build indices at official site of HISAT2.
    # https://ccb.jhu.edu/software/hisat2/manual.shtml
    input:
        # Required input.
        reference = ['reference/Homo_sapiens.GRCh38.genome.fasta.gz'],
        annotation = 'reference/gencode.v28.gtf.gz'
    output:
        'reference/hisat2_index/Homo_sapiens.GRCh38.genome.1.ht2'
    threads: 1
    params:
        # Optional parameters. Read through the comments carefully.
        # Provide appropriate option for your data,
        # or comment out the option if it is not needed.
        extra = '',
        # SNP file name
        snp = '',
        haplotype = '',
        splice_site = '',
        exon = '',
    log: 'logs/hisat/build/hisat-build.log'
    wrapper:
        'http://dohlee-bio.info:9193/hisat2/build'
