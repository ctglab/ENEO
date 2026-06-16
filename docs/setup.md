# Workflow Setup 

The pipeline uses the configuration file hosted in the `config` directory, named `config_main.yaml`. It's a YAML file divided into sections, that is actively read by Snakemake to obtain information about resources, parameters and paths. While we tried to make it as user-friendly as possible, enabling the automatic update during the downloading of the resources (refer to the [Setup resources](https://ctglab.github.io/ENEO/resources) section), some tweaks are required to do manually.

The configuration file is divided into sections, and here we provide a detailed explanation of each section and its content.

## General

The first three parameters are put on top of the configuration file, and they are the most important ones. They are:

- `execution_mode`: it's just a placeholder. Leave it with `full`.
- `OUTPUT_FOLDER`: the path to the output folder where the results will be saved. It is mandatory to provide the **absolute** path.
- `TEMP_DIR`: the path to the temporary folder where the intermediate files will be saved. It is mandatory to provide the **absolute** path.

```yaml
execution_mode: "full"
OUTPUT_FOLDER: /path/to/../ENEO_output/
TEMP_DIR: /path/to/../ENEO_temp/
```

!!! note
    Both the `OUTPUT_FOLDER` and `TEMP_DIR` paths must be absolute, and they're placed on top as they must be provided as mounting points in the `SINGULARITY_ARGS` parameters Snakemake.

## Datadirs

Each entry in this section defines the name of the output folders that will be created in the `OUTPUT_FOLDER` path. The name of the folder is the key, and the value is the path to the folder. The pipeline will create the folders if they don't exist, and will save the results in the corresponding folder.
<u> We suggest to keep the default values </u>, but you can change them if you want.

| Folder Name      | Description                                                                     |
|------------------|---------------------------------------------------------------------------------|
| `bams`           | intermediate BAM files during processing.                                       |
| `BQSR`           | recalibrated BAM files after base quality score recalibration.                  |
| `HLA_typing`     | HLA typing results.                                                             |
| `index_folder`   | STAR genome index.                                                              |
| `mapped_reads`   | STAR alignment output.                                                          |
| `peptides`       | pMHC binding affinity predictions and filtered epitope tables.                  |
| `qc_reports`     | aggregated quality control reports.                                             |
| `salmon_idx`     | Salmon transcript index.                                                        |
| `salmon_quant`   | per-sample Salmon quantification output.                                        |
| `expression`     | expression data derived from Salmon quantification.                             |
| `trimmed_reads`  | fastp-trimmed reads and rRNA-depleted reads (intermediate, removed after use).  |
| `trimming_report`| fastp HTML and JSON quality reports.                                            |
| `utils`          | workflow utilities.                                                             |
| `VCF`            | raw variant calls.                                                              |
| `VCF_out`        | processed and filtered VCF files.                                               |

## Parameters

This section contains the parameters used by the various steps of the pipeline. The majority of them are standard and should not be changed.

### Preprocessing

Reads are trimmed with `fastp` before rRNA depletion with `SortMeRNA`.

```yaml
params:
  fastp:
    threads: 6
    extra: "-q 20 -u 20 -l 50 -y 20 -x -g -3 -e 30 --detect_adapter_for_pe"
  sortmerna:
    threads: 8
```

| Parameter | Description |
|-----------|-------------|
| `fastp.extra` | Command-line flags passed directly to fastp. The defaults enforce a minimum base quality of 20 (`-q`), discard reads shorter than 50 bp (`-l`), and enable polyX tail trimming (`-x`). Refer to the [fastp documentation](https://github.com/OpenGene/fastp) for the full list of options. |
| `sortmerna.threads` | Number of threads allocated to the rRNA depletion step. |

### Variant calling

Variant calling is performed in parallel by Strelka2 and DeepVariant. Only variants called by both tools are retained.

```yaml
params:
  deepvariant:
    threads: 4
    extra: "split_skip_reads=true,channels=''"
```

| Parameter | Description |
|-----------|-------------|
| `deepvariant.threads` | Number of shards used by DeepVariant for parallel processing (`--num_shards`). |
| `deepvariant.extra` | Arguments forwarded to `make_examples`. `split_skip_reads=true` enables handling of split reads in RNA-seq data. |

### Neoantigen prediction

```yaml
params:
  pMHC:
    min_length: 8
    max_length: 12
    germProb: 0.5
```

| Parameter | Description |
|-----------|-------------|
| `min_length` | Minimum length of the mutated peptide submitted to NetMHCpan. |
| `max_length` | Maximum length of the mutated peptide submitted to NetMHCpan. |
| `germProb`   | Maximum germline probability allowed for a variant to generate a candidate peptide. Variants above this threshold are discarded. |

## Resources

This section contains the paths to the resources used by the pipeline. Most of the resources are downloaded and configured automatically by `setup/download_res.py`, which updates this section in place (refer to the [Setup resources](https://ctglab.github.io/ENEO/resources) section).

!!! note
    All paths must be **absolute**. It is generally preferred to keep all resources under the same root folder, as that folder must be mounted into the Singularity container via `--singularity-args`. If resources are spread across multiple locations, each one must be listed as a separate bind mount.

The following table summarises the resources expected by the pipeline and whether the setup script handles them automatically.

| Key | Description | Auto-downloaded |
|-----|-------------|-----------------|
| `genome` | GRCh38 reference genome (GIAB masked version). | yes |
| `transcriptome` | Gencode v47 transcript sequences (FASTA), used for Salmon indexing. | yes |
| `gtf` | Gencode v47 primary assembly annotation (GTF), used for STAR indexing. | yes |
| `dbsnps` | dbSNP population frequency VCF (ALFA release). | yes |
| `gsnps` | 1000 Genomes phase 1 high-confidence SNPs. | yes |
| `indel` | GATK known indels resource. | yes |
| `gnomad` | gnomAD allele frequency VCF (somatic hg38 subset). | yes |
| `PoN` | 1000 Genomes panel of normals for somatic filtering. | yes |
| `REDI` | REDI portal RNA-editing sites (BED format). | yes |
| `vep_cache` | VEP offline cache (v105, GRCh38). | yes |
| `sortmerna_db` | SortMeRNA default rRNA database (`smr_v4.3_default_db.fasta`). | yes |
| `deepvariant_rna_model` | DeepVariant RNA model directory (v1.4.0, inception_v3, RNA-seq standard). | yes |
| `giab_intervals` | GIAB hard-to-call regions used to filter variant calls (bundled). | no |
| `intervals_coding` | Exonic intervals for protein-coding genes used during variant calling (bundled). | no |


