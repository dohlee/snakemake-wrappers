from pathlib import Path
DATA_DIR = Path('YOUR_DATA_DIRECTORY')
RESULT_DIR = Path('YOUR_RESULT_DIRECTORY')

rule homer_motif_enrichment:
    input:
        bed = DATA_DIR / 'YOUR_INTERVAL_FILE_GOES_HERE',
        genome = 'YOUR_GENOME_FASTA_GOES_HERE',
    output:
        homer_results = directory(str(RESULT_DIR / '{sample}' / 'homerResults')),
        known_results = directory(str(RESULT_DIR / '{sample}' / 'knownResults')),
    params:
        # Mask repeats/lower case sequence.
        # Default: False
        mask = False,
        # Removes background positions overlappig with target positions.
        # Default: False
        bg = False,
        # Chop up large background regions to the avg size of target regions.
        # Default: False
        chopify = False,
        # Motif length.
        # NOTE: values greater than 12 may cause the program to run out of memory.
        # In these cases decrease the number of sequences analyzed (-N),
        # or try analyzing shorter sequence regions (i.e. -size 100).
        # Default: 8,10,12
        len_ = '8,10,12',
        # Fragment size to use for motif finding.
        # Options:
        # -size <#,#> : For example, -size -100,50 will get sequences from -100 to +50 relative from center.
        # -size given : uses the exact regions you give it.
        # Default: 200
        size = 200,
        # Numter of motifs to optimize.
        # Default: 25
        S = 25,
        # Global optimization, searches for strings with # mismatches.
        # Default: 2
        mis = 2,
        # Don't search reverse strand for motifs.
        # Default: False
        norevopp = False,
        # Don't search for de novo motif enrichment.
        # Default: False
        nomotif = False,
        # Output RNA motif logos and compare to RNA motif database, automatically sets -norevopp.
        # Default: False
        rna = False,
        # Known motif options/visualization
        # Check against motif collections.
        # Default: auto
        mset = False,
        # Just visualize de novo motifs, don't check similarity with known motifs.
        # Default: False
        basic = False,
        # Scale sequence logos by information content.
        # Default: Doesn't scale.
        bits = False,
        # Don't search for de novo vs. known motif similarity.
        # Default: False
        nocheck = False,
        # Known motifs to check against de novo motifs.
        # Default: False,
        mcheck = False,
        # [DANGEROUS] Allow adjustment of the degeneracy threshold for known motifs to improve p-value.
        # Default: False
        float_ = False,
        # Don't search for known motif enrichment.
        # Default: False
        noknown = False,
        # Known motifs to check for enrichment.
        # Default: False
        mknown = False,
        # Omit humor.
        # Default: False
        nofacts = False,
        # Use weblogo/seqlogo/ghostscript to generate logos, default uses SVG now.
        # Default: False
        seqlogo = False,
        # Sequence normalization options
        # Use GC% for sequence content normalization, now the default.
        # Default: True
        gc = True,
        # Use CpG% instead of GC% for sequence content normalization.
        # Default: False
        cpg = False,
        # No CG correction.
        # Default: False
        noweight = False,
        # Advanced options.
        # Use hypergeometric for p-values, binomial is default.
        # Default: False
        h = False,
        # Number of sequences to use for motif finding.
        # Default: max(50k, 2x input)
        N = False,
        # Use local background, # of equal size regions around peaks to use i.e. 2.
        # Default: False
        local = False,
        # Remove redundant sequences matching greater than # percent, i.e. redundant 0.5.
        # Default: False
        redundant = False,
        # Maximum percentage of N's in sequence to consider for motif finding.
        # Default: 0.7
        maxN = False,
        # Motifs to mask before motif finding.
        # -maskMotif <motif file 1> [motif file 2]
        # Default: False
        maskMotif = False,
        # Motifs to optimize or change length of
        # -opt <motif file 1> [motif file 2]...
        # Default: False
        opt = False,
        # Randomize target and background sequences labels.
        # Default: False
        rand = False,
        # Use file for target and background - first argument is list of peak ids for targets.
        # -ref <peak file>
        # Default: False
        ref = False,
        # Perform analysis of individual oligo enrichment.
        # Default: False
        oligo = False,
        # Dump fasta files for target and background sequences for use with other programs.
        # Default: False
        dumpFasta = False,
        # Force new background files to be created.
        # Default: False
        preparse = False,
        # Location to search for preparsed file and/or place new files.
        # Default: False
        preparsedDir = False,
        # Keep temporary files.
        # Default: False
        keepFiles = False,
        # Calculate empirical FDR for de novo discovery #=number of randomizations
        # Default: False
        fdr = False,
        # Use homer2 instead of original homer.
        # Default: True
        homer2 = True,
        # Length of lower-order oligos to normalize in background.
        # Default: 3
        nlen = 3,
        # Max normalization iterations.
        # Default: 160
        nmax = 160,
        # Weight sequences to neutral frequencies, i.e. 25%, 6.25%, etc.
        # Default: False
        neutral = False,
        # Lower-order oligo normalization for oligo table, use if -nlen isn't working well.
        # Default: False
        olen = False,
        # Maximum expected motif instance per bp in random sequence.
        # Default: 0.01
        e = 0.01,
        # Size in MB for statistics cache
        # Default: 500
        cache = 500,
        # Skip full masking after finding motifs, similar to original homer.
        # Default: False
        quickMask = False,
        # Stop looking for motifs when seed logp score gets above #.
        # Default: -10
        minlp = False,
    log:
        'logs/homer_motif_enrichment/{sample}.log'
    benchmark:
        repeat('benchmarks/homer_motif_enrichment/{sample}.tsv', 1)
    wrapper:
        'http://dohlee-bio.info:9193/homer/motif/enrichment'
