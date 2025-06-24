import os 

rule align:
    input:
        unpack(get_fastq),
        index=config["datadirs"]["index_folder"],
    output:
        bam=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["mapped_reads"],
                "{patient}_Aligned.out.bam"
            ),
    container:
        "docker://danilotat/eneo"
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
        runtime="960m",
        ncpus=4,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["align"],
            "{patient}.log"
        ),
    shell:
        """
        STAR --readFilesIn {input.r1} {input.r2} \
        --genomeDir {input.index} --runThreadN {threads} \
        --outFileNamePrefix {params.prefix} {params.extra}
        """


rule sortAlign:
    input:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["mapped_reads"],
            "{patient}_Aligned.out.bam"
        ),
    output:
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["mapped_reads"],
                "{patient}_Aligned.sortedByCoord.out.bam"
            ),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/samtools.yml" 
    threads: config["params"]["samtools"]["threads"]
    resources:
        mem="10G",
        runtime="120m",
        ncpus=2,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["sort_bam"],
            "{patient}.log"
        ),
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """


rule indexSortAligned:
    input:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["mapped_reads"],
            "{patient}_Aligned.sortedByCoord.out.bam"
        ),
    output:
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["mapped_reads"],
                "{patient}_Aligned.sortedByCoord.out.bam.bai"
            ),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    resources:
        mem="10G",
        runtime="60m",
        ncpus=2,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["sort_bam"],
            "{patient}.log"
        ),
    shell:
        """
        samtools index -@ {threads} {input}
        """
