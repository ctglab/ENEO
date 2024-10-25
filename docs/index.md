# Introduction

ENEO is a Snakemake workflow developed for the identification of cancer neoantigens using solely the tumor RNAseq, without requiring matched controls or additional sequencing experiments. You could read more from the preprint [here](https://www.biorxiv.org/content/10.1101/2024.08.08.607127v1).


## Quick Start

To start, clone the repo using 

```
git clone https://github.com/ctglab/ENEO.git
```

Then setup resources, as reported in the [dedicated section](https://ctglab.github.io/ENEO/resources).

The next step is to setup patients and their relative sequencing files. Edit the `units.csv` file to specify the path for fastq files of each patient.

    ```
    patient,fq1,fq2
    Pat_01,/path/to/Pat_01_1.fastq.gz,/path/to/Pat_01_2.fastq.gz
    Pat_02,/path/to/Pat_02_1.fastq.gz,/path/to/Pat_02_2.fastq.gz
    ```

Edit also the `patients.csv` file to add the list of patients to be processed

    ```
    patient
    Pat_01
    Pat_02
    ```

To execute the pipeline, be sure to have [snakemake](https://snakemake.readthedocs.io/en/stable/) and [singularity](https://docs.sylabs.io/guides/3.1/user-guide/index.html) installed. Then execute the pipeline using 

```
snakemake --use-singularity --use-conda --cores 4
```

If you spot any issue, please report in the [github issue section](https://github.com/ctglab/ENEO/issues)
