rule manta_tumor_normal:
    input:
        # Required input.
        tumor = '{tumor}.sorted.bam',
        tumor_index = '{tumor}.sorted.bam.bai',
        normal = '{normal}.sorted.bam',
        normal_index = '{normal}.sorted.bam.bai',
        reference = 'reference/Homo_sapiens_assembly38.fasta',
        reference_index = 'reference/Homo_sapiens_assembly38.fasta.fai'
    output:
        # Required output.
        candidate_small_indels = 'result/{tumor}/{tumor}_vs_{normal}.candidateSmallIndels.vcf.gz',
        candidate_sv = 'result/{tumor}/{tumor}_vs_{normal}.candidateSV.vcf.gz',
        diploid_sv = 'result/{tumor}/{tumor}_vs_{normal}.diploidSV.vcf.gz',
        somatic_sv = 'result/{tumor}/{tumor}_vs_{normal}.somaticSV.vcf.gz',
        candidate_small_indels_idx = 'result/{tumor}/{tumor}_vs_{normal}.candidateSmallIndels.vcf.gz.tbi',
        candidate_sv_idx = 'result/{tumor}/{tumor}_vs_{normal}.candidateSV.vcf.gz.tbi',
        diploid_sv_idx = 'result/{tumor}/{tumor}_vs_{normal}.diploidSV.vcf.gz.tbi',
        somatic_sv_idx = 'result/{tumor}/{tumor}_vs_{normal}.somaticSV.vcf.gz.tbi',
    params:
        extra = '',
        # Set options for WES input: turn off depth filters.
        # Default: False
        exome = False,
        # Set options for RNA-Seq input. Must specify exactly one bam input file.
        # Default: False
        rna = False,
        # Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand.
        # Default: False
        unstrandedRNA = False,
        # Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to
        # call. No VCF output will be provided outside of these regions. The full genome will be provided
        # outside of these regions. The full genome will still be used to estimate statistics from the input
        # (such as expected fragment size distribution). Only one BED file may be specified.
        # Default: call the entire genome.
        callRegions = False,
        # Pre-calculated alignment statistics file. Skips alignment stats calculation.
        # Default: None
        existingAlignStatsFile = False,
        # Use pre-calculated chromosome depths.
        # Default: None
        useExistingChromDepths = False,
        # Keep all temporary files (for workflow debugging).
        # Default: False
        retainTempFiles = False,
        # Generate a bam of supporting reads for all SVs.
        # Default: False
        generateEvidenceBam = False,
        # Output assembled contig sequences in VCF file.
        # Default: False
        outputContig = False,
        # Maximum sequence region size (in megabases) scanned by each task during SV Locus graph generation.
        # Default: 12
        scanSizeMb = 12,
    threads: 4
    log: 'logs/manta/tumor-normal/{tumor}_vs_{normal}.log'
    benchmark: 'benchmarks/manta/tumor-normal/{tumor}_vs_{normal}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/manta/tumor-normal'
