# Setup resources

ENEO heavily depends on public genetic databases for germline probability estimation plus other fairly common resources daily needed for bioinformatics. In order to make the download and configuration of resources for the workflow as smooth as possible, a python configuration script is available inside `setup/download_res.py`. As the workflow uses files with different annotations, downloading resources is not enough, as they may need a chromosome naming conversion: thus some common used tools like `bcftools` are required. To run the automatic setup, first create a `conda` environment using the `setup_env.yml` file located inside the `setup` folder.

```sh
conda env create -f setup_env.yml 
```

Then activate the conda environment, and launch the configuration script as follows, replacing the `resources` folder with .

```bash
python3 download_res.py ../config/config_main.yaml resources.json resources  
```

!!! note
    The genome assembly in use is the GRCh38. We're not planning to backport it to older assembly. 


This script will download different resources:
- genome, transcriptome and GTF from [Gencode](https://www.gencodegenes.org/human/)
- 1000G, ExAC from the GATK resource bundle
- dbSNPs ALFA 
- REDI portal
- VEP cache (v105)

## What if I already got some of them?

The configuration script works by controlling the existence of the files whose path is written inside the main configuration file `conf_main.yaml`, located in the `config` folder. If any of those files are already in your machine, just edit the configuration file adding the right *absolute* path. The script will check for its presence without re-downloading it.

!!! warning
    The annotation for files must be concordant throughout the pipeline, to avoid any error. Be sure to add resources manually only if you're sure of their concordance.


