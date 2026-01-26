# This is a common rule file intended to be used just for
# Github CI testing

import pandas as pd
import os
import glob
from pathlib import Path
from snakemake.utils import min_version

min_version("8.0.0")


configfile: "config/config.yaml"


configpath = "config/config.yaml"

patients = pd.read_csv("config/patients.csv")["patient"]
units = pd.read_csv("config/units.csv").set_index(["patient"], drop=False)
units = units.sort_index()
execution_mode = config["execution_mode"]
slurm_logdir = config["slurm_log_dir"]
logpath = Path(slurm_logdir)
logpath.mkdir(parents=True, exist_ok=True)
germProb_script = config["resources"]["germline_prob_script"]
bam_final_path = config["datadirs"]["BQSR"]
ref_fasta = config["resources"]["genome"]
ref_dict = ref_fasta.replace(".fa.gz", ".dict")
intervals_path = os.path.join(
    config["OUTPUT_FOLDER"] + config["datadirs"]["utils"], "interval-files"
)

num_workers = 20

READ = ["1", "2"]


wildcard_constraints:
    patient="|".join(patients),


def sample_from_patient(df, patient_list, condition):
    samples = []
    for x in patient_list:
        samples.append(
            df[(df.phenotype == condition) & (df.subject_id == x)].Sample.values[0]
        )
    return samples


def get_fastq(wildcards):
    """Return a dict where keys are read1/2 and values are list of files"""
    return {
        "r1": units.loc[wildcards.patient, "fq1"],
        "r2": units.loc[wildcards.patient, "fq2"],
    }


def memory_for_gatk(gatk_mem: int):
    """Quick workaround to return string parameter for gatk"""
    as_str = f'--java-options "-Xmx{gatk_mem}g"'
    return as_str


def get_bam(wildcards):
    """
    Return list of bam files for rule bam_readcount
    """
    return [f"{bam_final_path}/{wildcards.patient}_recal.bam"]


def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = "0" * num_zeros + str(i)
        ints.append(interval)
    return ints


def get_interval_files():
    ints = get_intervals()
    files = [i + "-scattered.interval_list" for i in ints]
    files = [os.path.join(intervals_path, f) for f in files]
    return files


def sample_from_patient(df, patient_list, condition):
    samples = []
    for x in patient_list:
        samples.append(
            df.loc[(df["phenotype"] == condition) & (df["subject_id"] == x)].values[0]
        )
    return samples


interval_files = get_interval_files()
