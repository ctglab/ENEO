import os 

rule align:
    input:
        r1=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["trimmed_reads"],
            "{patient}_1.fastq.gz"
        ),
        r2=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["trimmed_reads"],
            "{patient}_2.fastq.gz"
        ),
        index=config["datadirs"]["index_folder"],
    output:
        bam=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["mapped_reads"],
                "{patient}_Aligned.out.bam"
            )),
        star_log=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["mapped_reads"],
                "{patient}_Log.final.out"
            )),
    container:
        "docker://ctglabcnr/star"
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
        ncpus=1,
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
        "docker://ctglabcnr/eneo"
    conda:
        "../envs/samtools.yml" 
    params:
        threads=config["params"]["samtools"]["threads"]
    resources:
        mem="10G",
        runtime="120m",
        ncpus=1,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["sort_bam"],
            "{patient}.log"
        ),
    shell:
        """
        samtools sort -@ {params.threads} -o {output} {input}
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
        "docker://ctglabcnr/eneo"
    conda:
        "../envs/samtools.yml"
    params:
        threads=config["params"]["samtools"]["threads"]
    resources:
        mem="10G",
        runtime="60m",
        ncpus=1,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["sort_bam"],
            "{patient}.log"
        ),
    shell:
        """
        samtools index -@ {params.threads} {input}
        """
