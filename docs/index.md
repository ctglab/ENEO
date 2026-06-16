# Introduction

ENEO is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow developed for the identification of cancer neoantigens using solely the tumor RNAseq, without requiring matched controls or additional sequencing experiments. You could read more from the publication [here](https://academic.oup.com/nargab/article/7/3/lqaf026/8196479).

## Pipeline overview

Raw reads are first trimmed with [fastp](https://github.com/OpenGene/fastp) to remove adapters and low-quality bases, then passed through [SortMeRNA](https://github.com/sortmerna/sortmerna) to deplete ribosomal RNA reads before alignment. After STAR alignment and base quality score recalibration (BQSR), variant calling is performed in parallel by [Strelka2](https://github.com/Illumina/strelka) and [DeepVariant](https://github.com/google/deepvariant) using an RNA-specific model. Only variants concordantly called by both callers are retained and carried forward for annotation and neoantigen prediction.

## Quick Start

To execute the pipeline, it's required to have both [snakemake](https://snakemake.readthedocs.io/en/stable/) and [singularity](https://docs.sylabs.io/guides/3.1/user-guide/index.html) or [apptainer](https://apptainer.org/docs/user/latest/) installed. The easiest way to install both of them is using a dedicated conda environment. To create a new environment, use the following commands:

```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake apptainer
```

To start, clone the repo using 

```
git clone https://github.com/ctglab/ENEO.git
```

The next step is to setup resources, as reported in the [dedicated section](https://ctglab.github.io/ENEO/resources). 

Then update the configuration file in `config/config_main.yaml` following the instructions in the [dedicated section](https://ctglab.github.io/ENEO/setup).

The next step is to define the patients and their sequencing files. Edit `units.csv` to specify the **absolute** paths to the paired-end FASTQ files for each patient:

```
patient,fq1,fq2
Pat_01,/path/to/Pat_01_1.fastq.gz,/path/to/Pat_01_2.fastq.gz
Pat_02,/path/to/Pat_02_1.fastq.gz,/path/to/Pat_02_2.fastq.gz
```

Edit also `patients.csv` to list the patients to be processed. All entries must match the patient identifiers in `units.csv`.

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

## Reporting issue

If you encounter any error in the pipeline setup/execution, or you have any question regarding its usage, fill an [issue on github](https://github.com/ctglab/ENEO/issues) 
