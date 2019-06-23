__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2019, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Execute shell command.
shell(
    "("
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/HBR1.read1.fastq.gz -qO {snakemake.output[0]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/HBR1.read2.fastq.gz -qO {snakemake.output[1]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/HBR2.read1.fastq.gz -qO {snakemake.output[2]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/HBR2.read2.fastq.gz -qO {snakemake.output[3]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/HBR3.read1.fastq.gz -qO {snakemake.output[4]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/HBR3.read2.fastq.gz -qO {snakemake.output[5]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/UHR1.read1.fastq.gz -qO {snakemake.output[6]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/UHR1.read2.fastq.gz -qO {snakemake.output[7]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/UHR2.read1.fastq.gz -qO {snakemake.output[8]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/UHR2.read2.fastq.gz -qO {snakemake.output[9]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/UHR3.read1.fastq.gz -qO {snakemake.output[10]} && "
    "wget https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/test-data/rna-seq/pe/UHR3.read2.fastq.gz -qO {snakemake.output[11]}"
    ")"
)
