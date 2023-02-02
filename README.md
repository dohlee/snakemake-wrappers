# snakemake-wrappers

![pipe](https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/img/pipelines_cropped.jpeg)

Snakemake wrappers for bioinformatics research.

<a href="https://snakemake.bitbucket.io"><img src="https://img.shields.io/badge/snakemake-≥5.2.2-brightgreen.svg?style=flat-square" /></a>

## Available wrappers

<table style="border: 2px;">
  <tr>
    <td> Tool </td>
    <td> Description </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Differential expression</i> </td>
  </tr>
  <tr>
    <td> DESeq2 </td>
    <td> </td>
  </tr>
  <tr>
    <td> EBSeq </td>
    <td> </td>
  </tr>
  <tr>
    <td> edgeR </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Gene expression quantification</i> </td>
  </tr>
  <tr>
    <td> <a href="rsem">rsem</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="feature-counts">feature-counts</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="htseq">htseq</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="kallisto">kallisto</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="cufflinks">cufflinks</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="salmon">salmon</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Variant calling</i> </td>
  </tr>
  <tr>
    <td> <a href="freebayes">freebayes</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="somatic-sniper">somatic-sniper</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="varscan/somatic">varscan/somatic</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="strelka2">strelka2</a> </td>
    <td> </td>
  </tr>
  
  <tr style="background-color:#333333">
    <td colspan="2"> <i>Variant annotation</i> </td>
  </tr>
  <tr>
    <td> <a href="snpeff">snpeff</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Copy number analysis</i> </td>
  </tr>
  <tr>
    <td> <a href="varscan/copynumber">varscan/copynumber</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="varscan/copycaller">varscan/copycaller</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="dnacopy/cbs-segmentation">dnacopy/cbs-segmentation</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="cnvkit">cnvkit</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Alignment manipulation</i> </td>
  </tr>
  <tr>
    <td> <a href="sambamba">sambamba</a> </td>
    <td> A faster alternative of samtools. Wrappers for Sort, Index, Merge, Flagstat, Markdup are implemented for now.</td>
  </tr>
  <tr>
    <td> <a href="samtools">samtools</a> </td>
    <td> Manipulates Sequence Alignment/Map (SAM) format. Wrappers for Sort, Index, Faidx are implemented for now.</td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Sequencing data retrieval</i> </td>
  </tr>
  <tr>
    <td> <a href="sra-tools">sra-tools</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="prefetch">prefetch</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="parallel-fastq-dump">parallel-fastq-dump</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Preprocessing</i> </td>
  </tr>
  <tr>
    <td> <a href="fastqc">fastqc</a> </td>
    <td> The most popular quality-control tool for high-throughput sequencing data. </td>
  </tr>
  <tr>
    <td> <a href="trim-galore">trim-galore</a> </td>
    <td> Adapter trimming and quality control. Also provides RRBS-specific trimming option. </td>
  </tr>
  <tr>
    <td> <a href="fastp">fastp</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Sequencing read alignment</i> </td>
  </tr>
  <tr>
    <td> <a href="bwa">bwa</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="bowtie">bowtie</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="bowtie2">bowtie2</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Splice-aware sequencing read alignment</i> </td>
  </tr>
  <tr>
    <td> <a href="hisat2">hisat2</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="star">star</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Bisulfite-treated sequencing read alignment</i> </td>
  </tr>
  <tr>
    <td> <a href="bismark">bismark</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="bwa-meth">bwa-meth</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="bsmap">bsmap</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Bisulfite-seq postprocessing</i> </td>
  </tr>
  <tr>
    <td> <a href="methyldackel">methyldackel</a> </td>
    <td> </td>
  </tr>

  <tr style="background-color:#333333">
    <td colspan="2"> <i>Miscellaneous</i> </td>
  </tr>
  <tr>
    <td> <a href="subsample-fastq">subsample-fastq</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="wgsim">wgsim</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="tabix">tabix</a> </td>
    <td> </td>
  </tr>
  <tr>
    <td> <a href="igvtools">igvtools</a> </td>
    <td> </td>
  </tr>

</table>