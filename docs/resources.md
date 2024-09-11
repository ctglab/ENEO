# Setup resources

ENEO heavily depends on public genetic databases for germline probability estimation plus other fairly common resources daily needed for bioinformatics. In order to make the download and configuration of resources for the workflow as smooth as possible, a python configuration script is available inside `setup/download_res.py`. As the workflow depends on Ensembl annotation, downloading resources is not enough, as they may need a chromosome naming conversion: thus some common used tools like `bcftools` are required. To run the automatic setup, first create a `conda` environment using the `setup_env.yml` file located inside the `setup` folder.

```sh
conda env create -f setup_env.yml 
```

Then activate the conda environment, and launch the configuration script as follows, specifying the right `outfolder`.

```bash
python3 download_res.py --config ../config/config_main.yaml --json resources.json --outfolder resources  
```

## What if I already got some of them?

The configuration script works by controlling the existence of the files whose path is written inside the main configuration file `conf_main.yaml`, located in the `config` folder. If any of those files are already in your machine, just edit the configuration file adding the right *absolute* path. The script will check for its presence without re-downloading it.  
⚠️**NOTE**⚠️: As stated before, the annotation for all the files is the one provided by the Ensembl consortium. Files like the `dbSNPs` are often available with different chromosome naming: include them manually only if you're sure that the naming is right.   
