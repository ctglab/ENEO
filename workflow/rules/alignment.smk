import os 

rule align:
    input:
        unpack(get_fastq),
        index=config["datadirs"]["index_folder"],
    output:
        bam=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["mapped_reads"]
            + "/"
            + "{patient}_Aligned.out.bam"
        ),
    conda:
        "../envs/star.yml"
    params:
        index=lambda wc, input: input.index,
        prefix=lambda wc, output: os.path.join(os.path.dirname(os.path.abspath(output.bam)), f"{wc.patient}_"),
        extra="--sjdbGTFfile {} {}".format(
            config["resources"]["gtf"], config["params"]["STAR"]["extra"]
        ),
    threads: config["params"]["STAR"]["threads"]
    resources:
        mem="60G",
        time="16:00:00",
        ncpus=4,
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["align"]
        + "/"
        + "{patient}.log",
    shell:
        """
        STAR --readFilesIn {input.r1} {input.r2} \
        --genomeDir {input.index} --runThreadN {threads} \
        --outFileNamePrefix {params.prefix} {params.extra}
        """


rule sortAlign:
    input:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["mapped_reads"]
        + "/"
        + "{patient}_Aligned.out.bam",
    output:
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["mapped_reads"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.bam"
        ),
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    resources:
        mem="10G",
        time="2:00:00",
        ncpus=2,
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["sort_bam"]
        + "/"
        + "{patient}.log",
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """


rule indexSortAligned:
    input:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["mapped_reads"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.bam",
    output:
        temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["mapped_reads"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.bam.bai"
        ),
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    resources:
        mem="10G",
        time="1:00:00",
        ncpus=2,
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["sort_bam"]
        + "/"
        + "{patient}.log",
    shell:
        """
        samtools index -@ {threads} {input}
        """
