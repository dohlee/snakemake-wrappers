rule dnacopy_cbs_segmentation:
    input:
        copynumber_calls = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.copynumber.called',
    output:
        segments = 'result/{tumor_sample}/{tumor_sample}_vs_{normal_sample}.segment'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/dnacopy/cbs-segmentation'
