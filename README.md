<h1 align="center">snakemake-wrappers</h1>
<p align="center">Snakemake wrappers for bioinformatics.</p>
<p align="center">
  <a href="https://snakemake.bitbucket.io"><img src="https://img.shields.io/badge/snakemake-≥5.2.2-brightgreen.svg?style=flat-square" /></a>
  <a href="https://circleci.com/gh/dohlee/snakemake-wrappers"><img src="https://circleci.com/gh/dohlee/snakemake-wrappers.svg?style=svg" /></a>
  <a href="https://www.codacy.com/app/apap950419/snakemake-wrappers?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=dohlee/snakemake-wrappers&amp;utm_campaign=Badge_Grade"><img src="https://api.codacy.com/project/badge/Grade/a09065e963d9486dbeafecf7dbe44da0" /></a>
</p>

<h2 align="center">Available wrappers</h2>

- Routine analysis
  - DEG discovery
    - DESeq2
    - EBSeq
    - edgeR
- Gene expression quantification
  - [rsem](rsem)
  - [feature-counts](feature-counts)
  - [htseq](htseq)
  - [kallisto](kallisto)
  - [cufflinks](cufflinks)
  - [salmon](salmon)
- Variant calling
  - [freebayes](freebayes)
  - [somatic-sniper](somatic-sniper)
  - [varscan/somatic](varscan/somatic)
  - [strelka2](strelka2)
- Variant annotation
  - [snpeff](snpeff)
- Copy number analysis
  - [varscan/copynumber](varscan/copynumber)
  - [varscan/copycaller](varscan/copycaller)
  - [dnacopy/cbs-segmentation](dnacopy/cbs-segmentation)
  - [cnvkit](cnvkit)
- Alignment manipulation
  - [sambamba](sambamba)
  - [samtools](samtools)
- Data retrieval
  - [sra-tools](sra-tools)
    - [prefetch](prefetch)
  - [parallel-fastq-dump](parallel-fastq-dump)
- Preprocessing
  - [fastqc](fastqc)
  - [trim-galore](trim-galore)
  - [fastp](fastp)
- Aligner
  - [bwa](bwa)
  - [bowtie](bowtie)
  - [bowtie2](bowtie2)
- Splice-aware aligner
  - [hisat2](hisat2)
  - [star](star)
- Bisulfite-aligner
  - [bwa-meth](bwa-meth)
  - [bismark](bismark)
  - [bsmap](bsmap)
- Bisulfite-seq postprocessing
  - [methyldackel](methyldackel)
- Misc
  - [subsample-fastq](subsample-fastq)
  - [wgsim](wgsim)
  - [tabix](tabix)
