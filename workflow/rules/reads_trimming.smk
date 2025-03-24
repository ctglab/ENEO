import os

# Althought they may seem equivalent, explicit extra parameters could be defined which will benefit of
# PE or SE sequencing.
rule trimming_pe:
    input:
        sample=[
            os.path.join(config["resources"]["FASTQ"], "{patient}_1.fastq.gz"),
            os.path.join(config["resources"]["FASTQ"], "{patient}_2.fastq.gz"),
        ],
    output:
        trimmed=[
            os.path.join(config["datadirs"]["trimmed_reads"], "{patient}_1.fastq.gz"),
            os.path.join(config["datadirs"]["trimmed_reads"], "{patient}_2.fastq.gz"),
        ],
        html=os.path.join(config["datadirs"]["trimming_report"], "{patient}_fastp.html"),
        json=os.path.join(config["datadirs"]["trimming_report"], "{patient}_fastp.json"),
    params:
        extra=config["params"]["fastp"]["pe"],
    threads: config["params"]["thread"],
    container:
        "docker://danilotat/eneo"
    log:
        os.path.join(config["datadirs"]["logs"]["trimming"], "{patient}.log"),
    wrapper:
        "v1.0.0/bio/fastp"
