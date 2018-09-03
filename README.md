<h1 align="center">snakemake-wrappers</h1>
<p align="center">Snakemake wrappers for bioinformatics.</p>
<p align="center">
  <a href="https://circleci.com/gh/dohlee/snakemake-wrappers"><img src="https://circleci.com/gh/dohlee/snakemake-wrappers.svg?style=svg" /></a>
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
