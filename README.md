<h1 align="center">snakemake-wrappers</h1>
<p align="center">Snakemake wrappers for bioinformatics.</p>
<p align="center">
  <a href="https://snakemake.bitbucket.io"><img src="https://img.shields.io/badge/snakemake-≥5.2.2-brightgreen.svg?style=flat-square" /></a>
  <a href="https://circleci.com/gh/dohlee/snakemake-wrappers"><img src="https://circleci.com/gh/dohlee/snakemake-wrappers.svg?style=svg" /></a>
  <a href="https://www.codacy.com/app/apap950419/snakemake-wrappers?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=dohlee/snakemake-wrappers&amp;utm_campaign=Badge_Grade"><img src="https://api.codacy.com/project/badge/Grade/a09065e963d9486dbeafecf7dbe44da0" /></a>
</p>

<h2 align="center">Available wrappers</h2>

### Routine analysis
- DEG discovery
  - DESeq2
  - EBSeq
  - edgeR

### Variant calling
- freebayes

### File format handlers
- sambamba
  - index
  - sort
  - markdup
- samtools
  - index
  - sort

### Data retrieval
- sra-tools
  - prefetch
    - accession
- parallel-fastq-dump

### Preprocessing
- fastqc
- trim-galore
- fastp

### Aligner
- bismark
  - genome-preparation
- bsmap

### Gene expression quantification
- rsem
  - calculate-expression
    - pe
  - prepare-reference
- kallisto
  - index
  - quant
- salmon
  - index
  - quant

### Misc
- subsample-fastq
