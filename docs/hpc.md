# Run on HPC

ENEO was developed and tested in High Performances Computing (HPC) clusters with the SLURM workload manager. Snakemake deeply changed the job submissions and handling after the major update introduced with the version 8.0.0. Currently more than a single way exists for submitting jobs using Snakemake, but the most effective one seems to be using the `cluster-generic` plugin. 

If you're using Snakemake > 8.0.0, install the cluster-generic plugin using pip

```
pip install snakemake-executor-plugin-cluster-generic
```

Then inside the folder `worflow/profile` you'll find for each of the supported method (SLURM/SGE) two configuration files: one with the string `v8` in the name, used by Snakemake version >8.0.0, and a legacy `config.yaml`, for older versions. 

!!! tip

The following notes reported examples using the legacy config file. However, the relevant edits are the same! 

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
cluster:
  mkdir -p slurm-logs/{rule} &&
  sbatch
    --cpus-per-task={resources.ncpus}
    --mem={resources.mem}
    --time={resources.time}
    --job-name=smk-{rule}-{wildcards}
    --output=slurm-logs/{rule}/{rule}-{wildcards}-%j.out
    --partition=<partitionhere>
    --account=<accounthere>
```

This will create a folder called `slurm-logs` with a subfolder for each rule, where each patient will have a different log file. 

Then execute the pipeline with

```
snakemake --profile workflow/profile/slurm_profile
```

## SGE

!!! warning
  The support for SGE is still experimental. If you spot any issue, report it in the Github section

A config file for SGE is under `workflow/profile/sge_profile/config.yaml`. The overall scheme is the following

```yaml
cluster:
  mkdir -p sge-logs/{rule} &&
  qsub
    -pe smp {resources.ncpus}
    -l mem_free={resources.mem}
    -l h_rt={resources.time}
    -N smk-{rule}-{wildcards}
    -o sge-logs/{rule}/{rule}-{wildcards}-$JOB_ID.out
    -e sge-logs/{rule}/{rule}-{wildcards}-$JOB_ID.err
    -q all.q
```

The behavior is analogous to the SLURM one. 

To execute the pipeline, run

```
snakemake --profile workflow/profile/sge_profile
```


