rule macs2_bdgcmp:
    input:
        # Required input.
        treatment = 'YOUR_TREATMENT_BEDGRAPH',
        control = 'YOUR_CONTROL_BEDGRAPH',
    output:
        # Output suffixes automatically determines the '-m' parameters.
        # '_ppois.bdg': use '-m ppois', Poisson P-value
        # '_qpois.bdg': use '-m qpois', Poisson Q-value
        # '_subtract.bdg': use '-m subtract', subtraction from treatment
        # '_FE.bdg': use '-m FE', linear scale fold enrichment
        # '_logFE.bdg': use '-m logFE', log10 fold enrichment (need to set pseudocount)
        # '_logLR.bdg': use '-m logLR', log10 likelihood between ChIP-enriched model
        # and open chromatin model (need to set pseudocount)
        # '_slogLR.bdg': use '-m slogLR', symmetric log10 likelihood between two ChIP-enrichment
        # models, or maximum value between the two tracks.
        ['{treat}_vs_{control}_ppois.bdg', '{treat}_vs_{control}_FE.bdg']
    params:
        extra = '',
        # Scaling factor for treatment and control track. Keep it as
        # 1.0 or default in most cases. Set it ONLY while you have SPMR
        # output from MACS2 callpeak, and plan to calculate scores as MACS2
        # callpeak module. If you want to simulate 'callpeak' w/o '--to-large',
        # calculate effective smaller sample size after filtering redundant
        # reads in million (e.g., put 31.415926 if effective reads are
        # 31,415,926) and input if for '-S'; for 'callpeak --to-large',
        # calculate effective reads in larger sample.
        # Default: 1.0
        scaling_factor = 1.0,
        # The pseudocount used for calculating logLR, logFE or FE.
        # The count will be applied after normalization of sequencing
        # depth.
        # Default: False
        pseudocount = False,
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_bdgcmp/{treatment}_vs_{control}.tsv', 1)
    log: 'logs/macs2_callpeak/{treatment}_vs_{control}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/bdgcmp'
