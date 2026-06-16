# Run on HPC

ENEO was developed and tested in High Performances Computing (HPC) clusters with the SLURM workload manager. It's strongly suggested to use a recent version of Snakemake (>8.0.0) that works smootly with the slurm executor plugin. 

Install the SLURM executor plugin with

```
pip install snakemake-executor-plugin-slurm
```

Then inside the folder `worflow/profile/slurm` you'll find a configuration file named `config.yaml`, where you should add the details about your SLURM account and desired partition.  

## Singularity args

Two rules of the workflow (variant annotation and pMHC binding affinity estimation) depend on Singularity containers. It's key to ensure that all the relevant folders are readable/writable within each container. For this reason, multiple folders are required to be mounted, as Snakemake is *lazy* in assigning mountpoints. 

Populate the last entry of the config file, `singularity-args`, adding the absolute path for:

 - the resources directory
 - the temporary directory
 - the output directory
 - the workflow directory

Additionally, you had to set the TMPDIR environment variable to the temporary directory, to avoid writing permissions in the last step.


## SLURM

Insert the account and partition inside `workflow/profile/slurm_profile/config.yaml` and any other additional flags required for submitting jobs on the HPC platform in use. 

``` yaml
executor: slurm
default-resources:
    slurm_account: 
    slurm_partition: 
    mem_mb_per_cpu: 8000
    runtime: "60m"

use-apptainer: true
keep-going: true
rerun-incomplete: true
jobs: 100
singularity-args: '-B /../ENEO_res -B /../ENEO_temp -B /../ENEO_output -B /../ENEO/workflow --env TMPDIR=/../ENEO_temp'
```

For additional guidelines on how to compile the profile and how to execute it under the SLURM scheduler, refer to the [executor plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) 

```
snakemake --profile workflow/profile/slurm_profile
```

