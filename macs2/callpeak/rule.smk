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
        peak = Path(output_dir) /  '{sample}.narrowPeak',
        excel = Path(output_dir) /  '{sample}_peaks.xls',
        summits = Path(output_dir) /  '{sample}_summits.bed',
        model_script = Path(output_dir) /  '{sample}_model.R',
    params:
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        genome_size = 'hs',

        # Call broad peaks?
        broad = False,

        # Optional parameters. Omit if unneeded.

        # Random seed for data downsampling.
        seed = 0,

        # Output as bedGraph format? (-B/--bdg)
        bedGraph_out = True,

        # Cutoff for minimum FDR.
        q_value_cutoff = 0.01,

        # Use cutoff for minimum p-value.
        # p_value_cutoff = 1e-5
        
        # Extra options.
        extra = '',
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat("benchmarks/macs2_callpeak/{sample}.tsv", 1)
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
        peak = Path(output_dir) /  '{sample}.broadPeak',
        excel = Path(output_dir) /  '{sample}_peaks.xls',
        summits = Path(output_dir) /  '{sample}_summits.bed',
        model_script = Path(output_dir) /  '{sample}_model.R',
    params:
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        genome_size = 'hs',

        # Call broad peaks?
        broad = True,

        # Optional parameters. Omit if unneeded.

        # Random seed for data downsampling.
        seed = 0,

        # Output as bedGraph format? (-B/--bdg)
        bedGraph_out = True,

        # Cutoff for minimum FDR.
        q_value_cutoff = 0.01,

        # Use cutoff for minimum p-value.
        # p_value_cutoff = 1e-5
        
        # Extra options.
        extra = '',
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat("benchmarks/macs2_callpeak/{sample}.tsv", 1)
    log: 'logs/macs2_callpeak/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'


