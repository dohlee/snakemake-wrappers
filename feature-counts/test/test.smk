from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()

rule all:
    input: 'test.feature_counts.result'

rule fetch_annotation:
    input: FTP.remote('ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz')
    output: 'gencode.v28.gtf.gz'
    shell: 'mv {input} {output}'

rule fetch_bam:
    input: HTTP.remote('sgp1.digitaloceanspaces.com/dohlee-bioinfo/test-data/bam/test.bam')
    output: 'test.bam'
    shell: 'mv {input} {output}'

rule feature_counts:
    input:
        # Required input.
        alignment = ['{sample}.bam'],
        annotation = 'gencode.v28.gtf.gz'
    output:
        # Required output.
        '{sample}.feature_counts.result'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Format of provided annotation file. SAF or GTF.
        annotation_format = 'GTF',
        # GTF-specific options.
        feature_type = 'exon',  # Feature type in GTF file.
        attribute_type = 'gene_id',  # Attribute type in GTF file.
        # Minimum number of overlapping bases in a read that is required for read assignment.
        # [default: 1]
        min_overlap = 1,
    threads: 4
    wrapper:
        'feature-counts'
