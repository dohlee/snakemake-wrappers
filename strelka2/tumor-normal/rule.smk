rule strelka2_tumor_normal:
    input:
        # Required input.
        tumor = '{tumor_sample}.sorted.bam',
        tumor_index = '{tumor_sample}.sorted.bam.bai',
        normal = '{normal_sample}.sorted.bam',
        normal_index = '{normal_sample}.sorted.bam.bai',
        reference = 'reference/Homo_sapiens_assembly38.fasta',
        reference_index = 'reference/Homo_sapiens_assembly38.fasta.fai'
    output:
        # Required output.
        snv = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.snvs.strelka.vcf.gz',
        indel = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.indels.strelka.vcf.gz'
    params:
        extra = '',
        # Output a bed file describing somatic callable regions of the genome.
        # Default: False
        outputCallableRegions = False,
        # Specify a VCF of candidate indel alleles. These alleles are always evaluated but only reported
        # in the output when they are inferred to exist in the sample. The VCF must be tabix indexed. All
        # indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored. This
        # option may be specified more than once, multiple input VCFs will be merged.
        # Default: False
        indelCandidates = False,
        # Specify a VCF of candidate alleles. These alleles are always evaluated and reported even if they
        # are unlikely to exist in the sample. The VCF must be tabix indexed. All indel alleles must be
        # left-shifted/normalzed, any unnormalized allele will trigger runtime error. This option may be 
        # specified more than once, multiple input VCFs will be merged. Note that for any SNVs provided in
        # the VCF, the SNV site will be reported (and for gVCF, excluded from block compression), but the
        # specific SNV alleles are ignored.
        # Default: None
        forcedGT = False,
        # Set options for exome or other targeted input: note in particular that this flag turns off
        # high-depth filters.
        # Default: False
        exome = False,
        # Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call.
        # No VCF output will be provided outside of these regions. The full genome will still be used to 
        # estimabe statistics from the input (such as expected depth per chromosome). Only one BED file
        # may be specified.
        # Default: False (call the entire genome)
        callRegions = False,
        # Extended options:
        # These options are either unlikely to be reset after initial site configuration or only of interest
        # for workflow development/debugging. They will not be printed here if a default exists unless
        # --allHelp is specified.
        # Noise vcf file (submit argument multiple times for more than one file)
        # Default: False
        noiseVcf = False,
        # Maximum sequence region size (in megabases) scanned by each task during genome variant calling. 
        # Default: 12
        scanSizeMb = 12,
        # Limit the analysis to one or more genome region(s) for debugging purposes. If this arugment is 
        # provided multiple times the union of all specified regions will be analyzed. All regions must be
        # non-overlapping to get meaningful result.
        # Examples: '--region chr20' (whole chromosome), '--region chr2:100-2000 --region chr3:2500-3000'
        # (two regions).
        # If this option is specified (one or more times) together with the --callRegions BED file,
        # then all region arguments will be intersected with the callRegions BED track.
        # Default: False
        region = False,
        # Set variant calling task memory limit (in megabytes). It is not recommended to change the default
        # in most cases, but this might be required for a sample of unusual depth.
        # Default: False
        callMemMb = False,
        # Keep all temporary files (for workflow debugging).
        # Default: False
        retainTempFiles = False,
        # Disable empirical variant scoring (EVS).
        # Default: False
        disableEVS = False,
        # Report all empirical variant scoring features in VCF output.
        # Default: False
        reportEVSFeatures = False,
        # Provide a custom empirical scoring model file for SNVs.
        # Default: False,
        snvScoringModelFile = False,
        # Provide a custome empirical scoring model file for indels.
        # Default: False
        indelScoringModelFile = False,
    threads: 4
    log: 'logs/strelka2/tumor-normal/{tumor_sample}_vs_{normal_sample}.log'
    benchmark: 'benchmarks/strelka2/tumor-normal/{tumor_sample}_vs_{normal_sample}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/strelka2/tumor-normal'
