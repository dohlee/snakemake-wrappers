from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
REMOTE_BASE = 'https://sgp1.digitaloceanspaces.com/dohlee-bioinfo/'

rule cnvkit_batch:
    input:
        tumor_bam = ['{tumor}.bam'],
        normal_bam = ['{normal}.bam'],
        reference = 'reference.fasta',
        # Available targets:
        # SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.bed
        # SureSelect_Human_All_Exon_V4.hg19.regions.bed
        # SureSelect_Human_All_Exon_V4+UTRs.hg19.regions.bed
        # SureSelect_Human_All_Exon_V5.hg19.regions.bed
        # SureSelect_Human_All_Exon_V5.hg38.regions.bed
        # SureSelect_Human_All_Exon_V5+UTRs.hg19.regions.bed
        # SureSelect_Human_All_Exon_V5+UTRs.hg38.regions.bed
        # SureSelect_Human_All_Exon_V6.hg19.regions.bed
        # SureSelect_Human_All_Exon_V6.hg38.regions.bed
        # SureSelect_Human_All_Exon_V6+UTRs.hg19.regions.bed
        # SureSelect_Human_All_Exon_V6+UTRs.hg38.regions.bed
        # SureSelect_Human_All_Exon_V7.hg19.regions.bed
        # SureSelect_Human_All_Exon_V7.hg38.regions.bed
        targets = HTTP.remote(REMOTE_BASE + 'resource/wxs-targets/SeqCap_EZ_Human_Exome_Library_V3.hg19.regions.bed', keep_local=True),
        # Available accessible regions:
        # access-5k-mappable.hg19.bed
        access = HTTP.remote(REMOTE_BASE + 'resource/accesible-regions/access-5k-mappable.hg19.bed', keep_local=True),
    output:
        output_reference = '{tumor}_vs_{normal}.cnn',
        tumor_ratios = '{tumor}.cnr',
        tumor_segments = '{tumor}.cns',
        normal_ratios = '{tumor}.cnr',
        normal_segments = '{tumor}.cns',
    params:
        # Tip: Use `--scatter` or `--diagram` to generate plots.
        extra = '',
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/cnvkit/batch'
