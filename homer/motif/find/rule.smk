from pathlib import Path
RESULT_DIR = Path('YOUR_RESULT_DIRECTORY')

rule homer_motif_find:
    input:
        bed = 'YOUR_INTERVAL_FILE_GOES_HERE',
        genome = 'YOUR_GENOME_FASTA_GOES_HERE',
        motif = 'YOUR_MOTIF_FILE_GOES_HERE',
    output:
        homer_result = RESULT_DIR / '{sample}' / '{sample}'.homer_result
    params:
        extra = ''
    log:
        'logs/homer_motif_find/{sample}.log'
    benchmark:
        repeat('benchmarks/homer_motif_find/{sample}.tsv', 1)
    wrapper:
        'http://dohlee-bio.info:9193/homer/motif/find'
