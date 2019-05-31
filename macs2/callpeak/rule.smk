from pathlib import Path
from os.path import splitext
output_dir = Path('output_directory')

rule macs2_callpeak_narrow:
    input:
        # Required input.
        treatment = 'YOUR_TREATMENT_SAMPLE',
        # Optional input.
        control = 'YOUR_CONTROL_SAMPLE',
    output:
        peak = Path(output_dir) /  '{sample}_peaks.narrowPeak',
        excel = Path(output_dir) /  '{sample}_peaks.xls',
        summits = Path(output_dir) /  '{sample}_summits.bed',
        model_script = Path(output_dir) /  '{sample}_model.r',
    params:
        # Extra options.
        extra = '',
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = 'hs',
        # It controls the MACS behavior towards duplicate tags at the
        # exact same location -- the same coordination and the same
        # strand. The 'auto' option makes MACS calculate the maximum tags
        # at the exact same location based on binomial distribution using
        # 1e-5 as pvalue cutoff; and the 'all' option keeps every tags. If
        # an integer is given, at most this number of tags will be kept at
        # the same location. Note, if you've used samtools or picard to flag reads
        # as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them
        # although the reads may be decided by MACS2 as duplicate later.
        # The default is to keep one tag at the same location.
        # Default: 1
        keep_dup = 1,
        # Buffer size for incrementally increasing internal array size to store
        # reads alignment information. In most cases, you don't have to change this
        # parameter. However, if there are large number of chromosomes/contigs/scaffolds
        # in your alignment, it's recommended to specify a smaller buffer size in order
        # to decrease memory useage (but it will take longer time to read alignment files).
        # Minimum memory requested for reading an alignment file is about # of CHROMOSOME *
        # BUFFER_SIZE * 2 Bytes.
        # Default: 100000
        buffer_size = 100000,
        # Experiment name, which will be used to generate output file names.
        # Default: "NA"
        name = False,
        # Whether or not to save extended fragment pileup, and local lambda tracks (two files)
        # at every bp into a bedGraph file.
        # Default: False
        bdg = False,
        # Tells MACS to include trackline with bedGraph files. To include this trackline while
        # displaying bedGraph at UCSC genome browser, can show name and description of the
        # file as well. However my suggestion is to convert bedGraph to bigWig, then show the
        # smaller and faster binary bigWig file at UCSC genome browser, as well as downstream
        # analysis. Require -B to be set.
        # Default: Not include trackline.
        trackline = False,
        # If True, MACS will save signal per million reads for fragment pileup profiles.
        # Require -B to be set.
        # Default: False
        SPMR = False,
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = False,
        # Band width for picking regions to compute fragment size. This value is only used
        # while building the shifting model.
        # Default: 300
        bw = 300,
        # Select the regions within MFOLD range of high-confidence enrichment ratio against
        # background to build model. Fold-enrichment in regions must be lower than upper limit,
        # and higher than the lower limit. Use as "-m 10 30"
        # Default: 5 50
        mfold = '5 50',
        # Whether turn on the auto pair model process. If se, when MACS failed to build paired model,
        # it will use the noodel settings, the --exsize parameter to extend each tags towards 3'
        # direction. Not to use this automate fixation is a default behavior now.
        # Default: False
        fix_bimodal = False,
        # Wheter or not to build the shifting model. If True, MACS will not build model.
        # By default it means shifting size = 100, try to set extsize to change it.
        # Default: False
        nomodel = False,
        # The arbitrary shift in bp. Use discretion while setting it other than default value. 
        # When NOMODEL is set, MACS will user this value to move cutting ends (5') towards 5'->3'
        # direction then apply EXTSIZE to extend them to fragments. When this value is negative,
        # ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for
        # ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting
        # enriched cuttin gloci such as certain DNAseI-Seq datasets. Note, you can't set values other
        # than 0 if format is BAMPE or BEDPE for paired-end data.
        # Default: 0
        shift = 0,
        # The arbitrary extension size in bp. When nomodel is tru, MACS will use this value as
        # fragment size to extend each read towards 3' end, then pile them up.
        # It's exactly twice then umber of obsolete SHIFTSIZE. In previous language,
        # each read is moved 5'->3' direction to middle of fragment by 1/2 d, then extended to
        # both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3'
        # into a d size fragment. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
        # Default: 200
        extsize = 200,
        # Minium FDR (q-value) cutoff for peak detection. -q and -p are mutually exclusive.
        # Default; 0.05
        qvalue = 0.05,
        # Pvalue cutoff for peak detection. -q and -p are mutually exclusive. If pvalue cutoff is
        # set, qvalue will not be calculated and reported as -1 in the final .xls file.
        # Default: not set.
        pvalue = False,
        # When set, scale the small sample up to the bigger sample. By default, the bigger dataset
        # will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and
        # more specific results. Keep in mind that scaling down will bring down background noise more.
        # Default: False
        to_large = False,
        # When set, use a custome scaling ratio of ChIP/control (e.g. calculated using NCIS) for 
        # linear scaling.
        # Default: ignore
        ratio = False,
        # When set, random sampling method will scale down the bigger sample. By default, MACS
        # uses linear scaling.
        # Default: False
        down_sample = False,
        # Set the random seed while down sampling data. Must be a non-negative integer.
        # Default: 0
        seed = 0,
        # Optional directory to store temp files.
        # Default: tmp
        tempdir = False,
        # If True, MACS will use fixed background lambda as local lambda for every peak region.
        # Normally, MACS calculates a dynamic local lambda to reflect the local bias due to 
        # potential chromatin structure.
        nolambda = False,
        # The small nearby region in basepairs to calculate dynamic lambda. This is used to
        # cap ture the bias near the peak summit region. Invalid if there is no control data.
        # If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will
        # always perform a d-size local lambda calculation. The final local bias should be the
        # maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 1000
        slocal = 1000,
        # The large nearby region in basepairs to calculate dynamic lambda. This is used to 
        # capture the surround bias. If you set this to 0, MACS will skip llocal lambda calculation.
        # *Note* that MACS will always perform a d-size local lambda calculation. The final local
        # bias should be the maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 10000
        llocal = 10000,
        # If set, MACS will try to call broad peaks by linking nearby highly enriched regions.
        # The linking region is controlled by another cutoff through --linking-cutoff.
        # The maximum linking region length is 4 times of d from MACS.
        # Default: False
        broad = False,
        # Cutoff for broad region. This option is not available unless --broad is set.
        # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff.
        # Default: 0.1
        broad_cutoff = 0.1,
        # While set, MACS2 will analyze number or total length of peaks that can be called
        # by different p-value cutoff then output a summary table to help user decide
        # a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file.
        # Note, minlen and maxgap may affect the results.
        # WARNING: May take ~30 folds longer time to finish.
        # Default: False
        cutoff_analysis = False,
        # If set, MACS will use a more sophisticated signal processing approach to find subpeak
        # summits in each enriched peak region.
        # Default: False
        call_summits = False,
        # When set, the value will be used to filter out peaks with low fold-enrichment.
        # Note, MACS2 use 1.0 as pseudocount while calculating fold-enrichment.
        # Default: 1.0
        fe_cutoff = 1.0,
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_callpeak/{sample}.tsv', 1)
    log: 'logs/macs2_callpeak/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'

rule macs2_callpeak_broad:
    input:
        # Required input.
        treatment = 'YOUR_TREATMENT_SAMPLE',
        # Optional input.
        control = 'YOUR_CONTROL_SAMPLE',
    output:
        peak = Path(output_dir) /  '{sample}_peaks.broadPeak',
        excel = Path(output_dir) /  '{sample}_peaks.xls',
        summits = Path(output_dir) /  '{sample}_summits.bed',
        model_script = Path(output_dir) /  '{sample}_model.r',
    params:
        # Extra options.
        extra = '',
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = 'hs',
        # It controls the MACS behavior towards duplicate tags at the
        # exact same location -- the same coordination and the same
        # strand. The 'auto' option makes MACS calculate the maximum tags
        # at the exact same location based on binomial distribution using
        # 1e-5 as pvalue cutoff; and the 'all' option keeps every tags. If
        # an integer is given, at most this number of tags will be kept at
        # the same location. Note, if you've used samtools or picard to flag reads
        # as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them
        # although the reads may be decided by MACS2 as duplicate later.
        # The default is to keep one tag at the same location.
        # Default: 1
        keep_dup = 1,
        # Buffer size for incrementally increasing internal array size to store
        # reads alignment information. In most cases, you don't have to change this
        # parameter. However, if there are large number of chromosomes/contigs/scaffolds
        # in your alignment, it's recommended to specify a smaller buffer size in order
        # to decrease memory useage (but it will take longer time to read alignment files).
        # Minimum memory requested for reading an alignment file is about # of CHROMOSOME *
        # BUFFER_SIZE * 2 Bytes.
        # Default: 100000
        buffer_size = 100000,
        # Experiment name, which will be used to generate output file names.
        # Default: "NA"
        name = False,
        # Whether or not to save extended fragment pileup, and local lambda tracks (two files)
        # at every bp into a bedGraph file.
        # Default: False
        bdg = False,
        # Tells MACS to include trackline with bedGraph files. To include this trackline while
        # displaying bedGraph at UCSC genome browser, can show name and description of the
        # file as well. However my suggestion is to convert bedGraph to bigWig, then show the
        # smaller and faster binary bigWig file at UCSC genome browser, as well as downstream
        # analysis. Require -B to be set.
        # Default: Not include trackline.
        trackline = False,
        # If True, MACS will save signal per million reads for fragment pileup profiles.
        # Require -B to be set.
        # Default: False
        SPMR = False,
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = False,
        # Band width for picking regions to compute fragment size. This value is only used
        # while building the shifting model.
        # Default: 300
        bw = 300,
        # Select the regions within MFOLD range of high-confidence enrichment ratio against
        # background to build model. Fold-enrichment in regions must be lower than upper limit,
        # and higher than the lower limit. Use as "-m 10 30"
        # Default: 5 50
        mfold = '5 50',
        # Whether turn on the auto pair model process. If se, when MACS failed to build paired model,
        # it will use the noodel settings, the --exsize parameter to extend each tags towards 3'
        # direction. Not to use this automate fixation is a default behavior now.
        # Default: False
        fix_bimodal = False,
        # Wheter or not to build the shifting model. If True, MACS will not build model.
        # By default it means shifting size = 100, try to set extsize to change it.
        # Default: False
        nomodel = False,
        # The arbitrary shift in bp. Use discretion while setting it other than default value. 
        # When NOMODEL is set, MACS will user this value to move cutting ends (5') towards 5'->3'
        # direction then apply EXTSIZE to extend them to fragments. When this value is negative,
        # ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for
        # ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting
        # enriched cuttin gloci such as certain DNAseI-Seq datasets. Note, you can't set values other
        # than 0 if format is BAMPE or BEDPE for paired-end data.
        # Default: 0
        shift = 0,
        # The arbitrary extension size in bp. When nomodel is tru, MACS will use this value as
        # fragment size to extend each read towards 3' end, then pile them up.
        # It's exactly twice then umber of obsolete SHIFTSIZE. In previous language,
        # each read is moved 5'->3' direction to middle of fragment by 1/2 d, then extended to
        # both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3'
        # into a d size fragment. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
        # Default: 200
        extsize = 200,
        # Minium FDR (q-value) cutoff for peak detection. -q and -p are mutually exclusive.
        # Default; 0.05
        qvalue = 0.05,
        # Pvalue cutoff for peak detection. -q and -p are mutually exclusive. If pvalue cutoff is
        # set, qvalue will not be calculated and reported as -1 in the final .xls file.
        # Default: not set.
        pvalue = False,
        # When set, scale the small sample up to the bigger sample. By default, the bigger dataset
        # will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and
        # more specific results. Keep in mind that scaling down will bring down background noise more.
        # Default: False
        to_large = False,
        # When set, use a custome scaling ratio of ChIP/control (e.g. calculated using NCIS) for 
        # linear scaling.
        # Default: ignore
        ratio = False,
        # When set, random sampling method will scale down the bigger sample. By default, MACS
        # uses linear scaling.
        # Default: False
        down_sample = False,
        # Set the random seed while down sampling data. Must be a non-negative integer.
        # Default: 0
        seed = 0,
        # Optional directory to store temp files.
        # Default: tmp
        tempdir = False,
        # If True, MACS will use fixed background lambda as local lambda for every peak region.
        # Normally, MACS calculates a dynamic local lambda to reflect the local bias due to 
        # potential chromatin structure.
        nolambda = False,
        # The small nearby region in basepairs to calculate dynamic lambda. This is used to
        # cap ture the bias near the peak summit region. Invalid if there is no control data.
        # If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will
        # always perform a d-size local lambda calculation. The final local bias should be the
        # maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 1000
        slocal = 1000,
        # The large nearby region in basepairs to calculate dynamic lambda. This is used to 
        # capture the surround bias. If you set this to 0, MACS will skip llocal lambda calculation.
        # *Note* that MACS will always perform a d-size local lambda calculation. The final local
        # bias should be the maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 10000
        llocal = 10000,
        # If set, MACS will try to call broad peaks by linking nearby highly enriched regions.
        # The linking region is controlled by another cutoff through --linking-cutoff.
        # The maximum linking region length is 4 times of d from MACS.
        # Default: False
        broad = True,
        # Cutoff for broad region. This option is not available unless --broad is set.
        # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff.
        # Default: 0.1
        broad_cutoff = 0.1,
        # While set, MACS2 will analyze number or total length of peaks that can be called
        # by different p-value cutoff then output a summary table to help user decide
        # a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file.
        # Note, minlen and maxgap may affect the results.
        # WARNING: May take ~30 folds longer time to finish.
        # Default: False
        cutoff_analysis = False,
        # If set, MACS will use a more sophisticated signal processing approach to find subpeak
        # summits in each enriched peak region.
        # Default: False
        call_summits = False,
        # When set, the value will be used to filter out peaks with low fold-enrichment.
        # Note, MACS2 use 1.0 as pseudocount while calculating fold-enrichment.
        # Default: 1.0
        fe_cutoff = 1.0,
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_callpeak/{sample}.tsv', 1)
    log: 'logs/macs2_callpeak/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'
