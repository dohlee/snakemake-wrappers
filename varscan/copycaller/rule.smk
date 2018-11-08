rule varscan_copycaller:
    input:
        # Required input.
        raw_copynumber_calls = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.copynumber'
    output:
        # Required output.
        copynumber_calls = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.copynumber.called',
        # Optional output.
        homdel = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.copynumber.called.homdel',
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Minimum normal read depth at a position to make a call [20]
        min_coverage = 20,
        # Minimum tumor read depth at a position to make a non-homdel call [10]
        min_tumor_coverage = 10,
        # Maximum depth in tumor for candidate homozygous deletions [5]
        max_homdel_coverage = 5,
        # Lower bound for log ratio to call amplification [0.25]
        amp_threshold = 0.25,
        # Upper bound for log ratio to call deletion (provide as positive number) [0.25]
        del_threshold = 0.25,
        # Minimum size (in bases) for a region to be counted [10]
        min_region_size = 10,
        # Recenter data around an adjusted baseline > 0 [0]
        recenter_up = 0
        # Recenter data around an adjusted baseline < 0 [0]
        recenter_down = 0
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/varscan/copycaller'
