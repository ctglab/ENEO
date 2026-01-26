import os

# Determine output filenames based on execution mode
# In CI mode, fastp outputs directly to final names (skip rRNA removal)
# In full mode, fastp outputs to intermediate files for sortmerna processing
if config.get("execution_mode") == "CI":
    _trimmed_r1_suffix = "{patient}_1.fastq.gz"
    _trimmed_r2_suffix = "{patient}_2.fastq.gz"
else:
    _trimmed_r1_suffix = "{patient}_trimmed_1.fastq.gz"
    _trimmed_r2_suffix = "{patient}_trimmed_2.fastq.gz"


rule trimming:
    input:
        unpack(get_fastq),
    output:
        r1=temp(
            os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["trimmed_reads"],
            _trimmed_r1_suffix
        )),
        r2=temp(
            os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["trimmed_reads"],
            _trimmed_r2_suffix
        )),
        html=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["trimmed_reads"],
            "{patient}_fastp.html"
        ),
        json=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["trimmed_reads"],
            "{patient}_fastp.json"
        ),
    params:
        extra=config["params"]["fastp"]["extra"],
        threads=config["params"]["fastp"]["threads"],
    resources:
        mem="20G",
        runtime="240m",
        ncpus=1,
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/fastp.yml"
    log:
        os.path.join(
            config["datadirs"]["logs"]["trimming"],
            "{patient}.log"),
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        -h {output.html} -j {output.json} \
        -w {params.threads} \
        {params.extra}
        """


# rRNA removal step - only included in full mode
if config.get("execution_mode") != "CI":
    rule remove_rrna:
        input:
            r1=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["trimmed_reads"],
                "{patient}_trimmed_1.fastq.gz"
            ),
            r2=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["trimmed_reads"],
                "{patient}_trimmed_2.fastq.gz"
            ),
            rrna_db=config["resources"]["sortmerna_db"],
        output:
            r1=temp(
                os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["trimmed_reads"],
                "{patient}_1.fastq.gz"
            )),
            r2=temp(
                os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["trimmed_reads"],
                "{patient}_2.fastq.gz"
            )),
            stats=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["trimmed_reads"],
                "{patient}_sortmerna.log"
            ),
        params:
            workdir=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["trimmed_reads"],
                "{patient}_sortmerna"
            ),
            out_prefix=lambda wc, output: os.path.join(
                os.path.dirname(output.r1), wc.patient
            ),
        threads: config["params"]["sortmerna"]["threads"]
        resources:
            mem="32G",
            runtime="240m",
            ncpus=1,
        container:
            "docker://danilotat/sortmerna:latest"
        conda:
            "../envs/sortmerna.yml"
        log:
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["logs"]["trimming"],
                "{patient}_sortmerna.log"
            ),
        shell:
            """
            sortmerna \
                --ref {input.rrna_db} \
                --reads {input.r1} \
                --reads {input.r2} \
                --workdir {params.workdir} \
                --aligned {params.workdir}/rrna \
                --other {params.out_prefix} \
                --paired_in \
                --fastx \
                --threads {threads} \
                --out2 2>&1 | tee {output.stats}
            mv {params.out_prefix}_fwd.fq.gz {output.r1}
            mv {params.out_prefix}_rev.fq.gz {output.r2}
            rm -rf {params.workdir}
            """
