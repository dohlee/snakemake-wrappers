rule bissnp_bisulfitegenotyper:
    input:
        # Input SAM or BAM file(s) for genotyping.
        bam = '{sample}.bam',
        # dbsnp VCF file (which should sort as the same order of chromosome in BAM file header), provide prior SNP information from dbSNP database. The chromosome contains order of dbSNP VCF file should be the same as your reference genome file and input BAM file's header.
        dbsnp = '/data/project/dohoon/resource/dbsnp151.vcf.gz',
        # Reference sequence fasta file (which should have the same chromosome's order as input BAM file's header). 
        reference = '/data/project/dohoon/reference/hg38/hg38.fa',
        # (Optional), A list of genomic intervals over which to operate. Can be explicitly specified on the commandline like '-L chr11:7000000-7100000', '-L chr1' or in a bed file.
        # intervals = '',
    output:
        # Output VCF file, when output  mode is DEFAULT_FOR_TCGA, this option is used to output all CpG sites. While -vnf2 option is used to output all sites. And -vnf2 option is required at that time. Required.
        vcf_file_name1 = '',
        # Output VCF file, when output mode is DEFAULT_FOR_TCGA, this option is required and used to store all SNP sites. In the other output mode, it is not required. Optional.
        vcf_file_name2 = '',
    params:
        extra = '',
        # The minimum phred-scaled threshold for genotype calling that is confident (which is marked PASS in VCF file)
        # Default: 20
        stand_call_conf = 20,
        # The minimum phred-scaled threshold for genotype calling that could be emitted.
        # Default: 0
        stand_emit_conf = 0,
        # Minimum read mapping quality required to consider a read for calling.
        # Default: 30
        min_mapping_quality_score = 30,
        # Minimum base mapping quality required to consider a base for calling.
        # Default: 17
        min_base_quality_score = 17,
        # What kind of output file do you want.
        # Options (Default: DEFAULT_FOR_TCGA)
        # EMIT_ALL_SITES: emit all of callable sites into vcf file.
        # EMIT_ALL_CONFIDENT_SITES: emit all of sites above emit confident threshold into vcf file.
        # EMIT_VARIANTS_ONLY: emit all of SNP sites above emit threshold into vcf file.
        # EMIT_ALL_CPG: emit all of CpG sites above emit threshold into vcf file.
        # EMIT_ALL_CYTOSINES: emit all of Cytosine sites above emit threshold into vcf file.
        # EMIT_HET_SNPS_ONLY: emit all of heterozygous SNP sites above emit threshold into vcf file.
        # EMIT_VARIANT_AND_CYTOSINES: emit all of Cytosine sites above emit threshold into vcf1 file, and all of SNP sites above emit threshold into vcf2 file.
        # DEFAULT_FOR_TCGA: emit all of CpG sites above emit threshold into vcf1 file, and all of SNP sites above emit threshold into vcf2 file.
        # 
        # Default: DEFAULT_FOR_TCGA
        output_modes = 'DEFAULT_FOR_TCGA',
    threads: 1
    log: 'logs/bissnp/bisulfitegenotyper/{sample}.log'
    benchmark: 'benchmarks/bissnp/bisulfitegenotyper/{sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/bissnp/bisulfitegenotyper'
