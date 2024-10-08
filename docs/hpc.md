# Run on HPC

ENEO was developed and tested in High Performances Computing (HPC) clusters with the SLURM workload manager. Even if Snakemake introduced plugins in version `>8.0`, still the preferred way to launch the workload in SLURM is using a defined profile.

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

Then execute the pipeline by just 

```
snakemake --profile workflow/profile/slurm_profile
```