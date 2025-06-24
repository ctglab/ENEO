import pandas as pd
import os
import glob
from snakemake.utils import min_version

min_version("5.9.1")


configfile: "config/config_main.yaml"

# Set execution mode
if config["execution_mode"].lower() == "reduced":
    execution_mode = "reduced"
    print(f"The pipeline will be executed in the {execution_mode} mode")
elif config["execution_mode"].lower() == "full":
    execution_mode = "full"
    print(f"The pipeline will be executed in the {execution_mode} mode")
else:
    raise ValueError(f"The parameter 'execution_mode' must be one of ['bam','fastq'].")

# Load patient info.
# Note that this dataframe is accessed every runtime to determine the wildcards used
# at each step of the analysis while needed.

configpath = "config/config_main.yaml"
patients = pd.read_csv("patients.csv")["patient"]
units = pd.read_csv("units_fixed.csv").set_index(["patient"], drop=False)
units = units.sort_index()

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


def keep_paired_samples(patient_list):
    counts = Counter(patient_list)
    # keep only patients with matched tumor/normal
    to_keep = [patient for patient in patient_list if counts[patient] > 1]
    # escamoutage to keep insertion order. Works only on Python > 3.6 !
    return list(dict.fromkeys(to_keep).keys())


def sample_from_patient(df, patient_list, condition):
    samples = []
    for x in patient_list:
        samples.append(
            df[(df.phenotype == condition) & (df.subject_id == x)].Sample.values[0]
        )
    return samples

def get_bam(wildcards):
    return {'bam': units.loc[wildcards.patient, 'bam']}

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
