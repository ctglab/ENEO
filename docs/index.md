# Introduction

ENEO is a Snakemake workflow developed for the identification of cancer neoantigens using solely the tumor RNAseq, without requiring matched controls or additional sequencing experiments. You could read more from the preprint [here](https://www.biorxiv.org/content/10.1101/2024.08.08.607127v1).


## Quick Start

To start, clone the repo using 

```
git clone https://github.com/ctglab/ENEO.git
```

To limit the size of the repository, test files are not provided directly while cloning but downloaded on-fly for CI testing. For proceeding using a *real* sample you could download tumor RNA-seq data from a patient of the NCI surgery branch from Steven Rosenberg group

``` bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR969/008/SRR9697628/SRR9697628_1.fastq.gz -O test_data/SRR9697628_1.fastq.gz && \
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR969/008/SRR9697628/SRR9697628_2.fastq.gz -O test_data/SRR9697628_2.fastq.gz
```

Then execute the pipeline using 

```
snakemake --use-singularity --use-conda --cores 4
```
If you spot any issue, please report in the github issue section https://github.com/ctglab/ENEO/issues
