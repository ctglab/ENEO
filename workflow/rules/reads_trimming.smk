# Althought they may seem equivalent, explicit extra parameters could be defined which will benefit of
# PE or SE sequencing.
rule trimming_pe:
    input:
        sample=[
            config["resources"]["FASTQ"] + "/" + "{patient}_1.fastq.gz",
            config["resources"]["FASTQ"] + "/" + "{patient}_2.fastq.gz",
        ],
    output:
        trimmed=[
            config["datadirs"]["trimmed_reads"] + "/" + "{patient}_1.fastq.gz",
            config["datadirs"]["trimmed_reads"] + "/" + "{patient}_2.fastq.gz",
        ],
        html=config["datadirs"]["trimming_report"] + "/" + "{patient}_fastp.html",
        json=config["datadirs"]["trimming_report"] + "/" + "{patient}_fastp.json",
    params:
        extra=config["params"]["fastp"]["pe"],
    threads: config["params"]["thread"]
    conda:
        "../envs/fastp.yml"
    log:
        config["datadirs"]["logs"]["trimming"] + "/" + "{patient}.log",
    wrapper:
        "v1.0.0/bio/fastp"
