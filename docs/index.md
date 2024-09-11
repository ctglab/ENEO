# Introduction

ENEO is a Snakemake workflow developed for the identification of cancer neoantigens using solely the tumor RNAseq, without requiring matched controls or additional sequencing experiments. You could read more from the preprint [here](https://www.biorxiv.org/content/10.1101/2024.08.08.607127v1).


## Quick Start

To start, clone the repo using 

```
git clone https://github.com/ctglab/ENEO.git
```

To execute the pipeline, be sure to have [snakemake](https://snakemake.readthedocs.io/en/stable/) and [singularity](https://docs.sylabs.io/guides/3.1/user-guide/index.html) installed. Then execute the pipeline using 

```
snakemake --use-singularity --use-conda --cores 4
```

If you spot any issue, please report in the [github issue section](https://github.com/ctglab/ENEO/issues)
