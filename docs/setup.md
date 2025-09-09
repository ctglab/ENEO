# Workflow Setup 

The pipeline uses the configuration file hosted in the `config` directory, named `config_main.yaml`. It's a YAML file divided into sections, that is actively read by Snakemake to obtain informations about resources, parameters and paths. While we tried to make it as user-friendly as possible, enabling the automatic update during the downloading of the resources (refer to the [Setup resources](https://ctglab.github.io/ENEO/resources) section), some tweaks are required to do manually.

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

| Folder Name    | Description                                                      |
|----------------|------------------------------------------------------------------|
| `bams`         | the name of the folder used during .BAM files processing.        |
| `BQSR`         | the name of the folder where the recalibrated .BAM files will be saved. |
| `HLA_typing`   | the name of the folder where the HLA typing results will be saved. |
| `index_folder` | the name of the folder where the STAR index will be saved.       |
| `mapped_reads` | the name of the folder where the results of STAR will be saved.  |
| `peptides`     | the name of the folder where the pMHC prediction will be saved.  |
| `utils`        | the name of the folder used for workflow utilities.              |
| `VCF`          | the name of the folder where the RAW .VCF files will be saved.   |
| `VCF_out`      | the name of the folder where the processed .VCF files will be saved. |

## Parameters

This section contains the parameters used by the various steps of the pipeline. The majority of the parameters are standard and should not be changed. We suggest anyway users to focus on the `pMHC` section of the parameters, which reports

```yaml
  pMHC:
    threads: 4
    netmhcpan_launcher_script: workflow/scripts/netmhcpan_launcher.py
    calibration_frame: workflow/supplementary_res/optimal_percentile_netmhcpan.csv
    hla_ligand_atlas: workflow/supplementary_res/HLA_ligand_atlas.tsv.gz
    filter_peptides_script: workflow/scripts/filter_peptides.py
    min_length: 8
    max_length: 12
    germProb: 0.5
```

The last three parameters defines the constrains used for the mutated peptide generation and the subsequent 

| Parameter Name | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `min_length`   | Minimum length of the peptide.                                         |
| `max_length`   | Maximum length of the peptide.                                         |
| `germProb`     | Threshold for the germline probability of the variant generating the mutated peptide.    |

## Resources

This section contains the paths to the resources used by the pipeline. Most of the resources could be downloaded automatically using the script provided in the folder `setup` named `download_res.py`. The path for the downloaded files will be updated automatically (refer to the [Setup resources](https://ctglab.github.io/ENEO/resources) section).

!!! note
    All the paths must be **absolute**. Note that it's generally preferred to keep all the resources in the same folder, but this is not mandatory. The reason is due to the fact that these files must be accessible by the Singularity container, and the paths must be provided as mounting points in the `SINGULARITY_ARGS` parameter when executing the pipeline. If all the files are in the same folder, you can provide the path to the folder, and the pipeline will be able to access all the files.


