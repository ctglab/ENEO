# Introduction

ENEO is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow developed for the identification of cancer neoantigens using solely the tumor RNAseq, without requiring matched controls or additional sequencing experiments. You could read more from the preprint [here](https://www.biorxiv.org/content/10.1101/2024.08.08.607127v1).


## Quick Start

To execute the pipeline, it's required to have both [snakemake](https://snakemake.readthedocs.io/en/stable/) and [singularity](https://docs.sylabs.io/guides/3.1/user-guide/index.html) installed. The easiest way to install both of them is using a dedicated conda environment. To create a new environment, use the following commands:

```
conda create -n eneo -c bioconda snakemake=8.6.0 
```

To start, clone the repo using 

```
git clone https://github.com/ctglab/ENEO.git
```

The next step is to setup resources, as reported in the [dedicated section](https://ctglab.github.io/ENEO/resources). 

Then update the configuration file in `config/config.yaml` following the instructions in the [dedicated section](https://ctglab.github.io/ENEO/setup).

The next step is to setup patients and their relative sequencing files, defined in the files `units.csv` and `patients.csv`. The pipeline is designed to be executed both from `FASTQ` or already aligned `BAM` files: this different behavior is specified using the `execution_mode` in the main configuration file. Edit then the `units.csv` file to specify the **absolute** path for sequencing files of each patient, accordingly:

- if you're executing the pipeline in `full` mode, you had to provide paths to the paired end fastq files, as detailed.

    ```
    patient,fq1,fq2
    Pat_01,/path/to/Pat_01_1.fastq.gz,/path/to/Pat_01_2.fastq.gz
    Pat_02,/path/to/Pat_02_1.fastq.gz,/path/to/Pat_02_2.fastq.gz
    ```

- alternatively, if you're executing the pipeline in `reduced` mode, you had to provide paths to the aligned bam files, as detailed.

    ```
    patient,bam
    Pat_01,/path/to/Pat_01.bam
    Pat_02,/path/to/Pat_02.bam
    ```
!!! note
    While running with BAM files save some time in the analysis, here the alignment is performed in Two-Pass mode to generate sample specific junctions. Thus it's suggested to run the pipeline in **full** mode. 

Edit also the `patients.csv` file to add the list of patients to be processed. All the listed patients must match the patients in the `units.csv` file.

    ```
    patient
    Pat_01
    Pat_02
    ```

## Executing the pipeline

If you're running the pipeline using an HPC cluster, we provided detailed execution profiles for SLURM schedulers under the `workflow/profile` directory and the relative instruction in the [dedicated section](https://ctglab.github.io/ENEO/hpc).

Alternatively, you can run the pipeline locally, using the following command:

```sh
snakemake \
--use-singularity \
--cores 8 \
--singularity-args "-B /path/to/ENEO -B /path/to/ENEO_output -B /path/to/ENEO_temp -B /path/to/eneo_resources --env TMPDIR=/path/to/ENEO_temp"
```
