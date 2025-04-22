rule BQSR_1:
    input:
        bam=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["bams"],
            "{patient}_split.out.bam"
        ),
        GSNPs=config["resources"]["gsnps"],
        indel=config["resources"]["indel"],
        DbSNP=config["resources"]["dbsnps"],
        fasta=config["resources"]["genome"],
    output:
        recall=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["bams"],
            "{patient}_recal.table"
        ),
    resources:
        time="6:00:00",
        ncpus=4,
        mem="32G",
    threads: config["params"]["BQSR"]["threads"]
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/gatk.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["base_recalibration"],
            "{patient}.log"
        ),
    shell:
        """
        gatk BaseRecalibrator \
        -I {input.bam} \
        -R {input.fasta} \
        --known-sites  {input.GSNPs} \
        --known-sites  {input.indel}  \
        -O {output.recall}
        """


rule applyBQSR:
    input:
        bam=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["bams"],
            "{patient}_split.out.bam"
        ),
        fasta=config["resources"]["genome"],
        recall=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["bams"],
            "{patient}_recal.table"
        ),
    output:
        rbam=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["BQSR"],
                "{patient}_recal.bam"
                ),
        )
    threads: config["params"]["BQSR"]["threads"]
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/gatk.yml"
    resources:
        time="6:00:00",
        ncpus=4,
        mem="32G",
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["base_recalibration"],
            "{patient}.log"
        ),
    shell:
        """
        gatk ApplyBQSR \
        -I {input.bam}  \
        -R {input.fasta} \
        --bqsr-recal-file {input.recall} \
        -O {output.rbam}
        """

rule compressBam:
    input:
        bam=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.bam"
        ),
        reference=config["resources"]["genome"],
    output:
        cram=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.cram"
        ),
        index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.cram.crai"
        ),
    threads: config["params"]["samtools"]["threads"],
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/samtools.yml"
    resources:
        time="2:00:00",
        ncpus=4,
        mem="32G",
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["base_recalibration"],
            "{patient}.log"
        ),
    shell:
        """
        samtools view -T {input.reference} -C -o {output.cram} {input.bam}
        samtools index {output.cram}
        """