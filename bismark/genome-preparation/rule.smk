rule bismark_genome_preparation:
    input:
        directory('reference/')
    output:
        directory('reference/Bisulfite_Genome')
    threads: 1
    wrapper:
        'http://dohlee-bio.info:9193/bismark/genome-preparation'
