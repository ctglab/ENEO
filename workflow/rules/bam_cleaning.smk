# these rules come from https://github.com/khandaud15/RNA-Seq-Variant-Calling/blob/master/Snakefile


rule AddGrp:
    input:
        bam=config["OUTPUT_FOLDER"]
        + config["datadirs"]["mapped_reads"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.bam",
    output:
        rg=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["bams"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.rg.bam"
        ),
    conda:
        "../envs/gatk.yml"
    params:
        RGPU="{patient}",
        RGSM="{patient}",
        tmp_dir=config["TEMP_DIR"],
        extra="--RGLB rg1 --RGPL illumina",
    resources:
        mem="32G",
        time="4:00:00",
        ncpus=4,
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["bam_cleaning"]
        + "/"
        + "{patient}.log",
    shell:
        """
        gatk AddOrReplaceReadGroups  -I {input.bam} -O {output.rg} {params.extra} --TMP_DIR {params.tmp_dir} --RGPU {params.RGPU} --RGSM {params.RGSM}
        """


rule bed_to_intervals:
    input:
        bed=config["resources"]["intervals_coding"],
        fasta_dict=ref_dict,
    output:
        intervals=config["OUTPUT_FOLDER"]
        + config["datadirs"]["utils"]
        + "/"
        + "coding.interval_list",
    conda:
        "../envs/gatk.yml"
    resources:
        time="1:00:00",
        ncpus=2,
        mem="8G",
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["intervals"]
        + "/"
        + "interval.log",
    shell:
        """
        gatk BedToIntervalList -I {input.bed} -SD {input.fasta_dict} -O {output.intervals}
        """

rule mark_duplicates:
    input:
        bam=config["OUTPUT_FOLDER"]
        + config["datadirs"]["bams"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.rg.bam",
    output:
        bam=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["bams"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.md.bam"
        ),
        metrics=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["bams"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.metrics.txt"
        ),
    params:
        hard_ram=config["params"]["gatk"]["RAM"],
        temporary_dir=config["TEMP_DIR"],
    conda:
        "../envs/gatk.yml"
    threads: config["params"]["MarkDuplicates"]["threads"]
    resources:
        mem="32G",
        time="4:00:00",
        ncpus=4,
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["bam_cleaning"]
        + "/"
        + "{patient}.log",
    shell:
        """
        gatk --java-options "-Xmx{params.hard_ram}g -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
        MarkDuplicates -I {input.bam} -O {output.bam} \
        -M {output.metrics} --ASSUME_SORT_ORDER coordinate --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
        --TMP_DIR {params.temporary_dir} 
        """


rule sort_bam_gatk:
    input:
        bam=config["OUTPUT_FOLDER"]
        + config["datadirs"]["bams"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.md.bam",
    output:
        bam_out=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["bams"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.md.sorted.bam"
        ),
    conda:
        "../envs/samtools.yml"
    resources:
        time="2:00:00",
        ncpus=2,
        mem="8G",
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["bam_cleaning"]
        + "/"
        + "{patient}.log",
    shell:
        """
        samtools sort {input.bam} -o {output.bam_out} 
        """


rule samtools_index:
    input:
        bam=config["OUTPUT_FOLDER"]
        + config["datadirs"]["bams"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.md.sorted.bam",
    output:
        bai=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["bams"]
            + "/"
            + "{patient}_Aligned.sortedByCoord.out.md.sorted.bam.bai"
        ),
    conda:
        "../envs/samtools.yml"
    resources:
        time="1:00:00",
        ncpus=2,
        mem="8G",
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["bam_cleaning"]
        + "/"
        + "{patient}.log",
    shell:
        """
        samtools index {input.bam} {output.bai} 
        """


rule SplitNCigarReads:
    input:
        bai=config["OUTPUT_FOLDER"]
        + config["datadirs"]["bams"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.md.sorted.bam.bai",
        bam=config["OUTPUT_FOLDER"]
        + config["datadirs"]["bams"]
        + "/"
        + "{patient}_Aligned.sortedByCoord.out.md.sorted.bam",
        intervals=config["OUTPUT_FOLDER"]
        + config["datadirs"]["utils"]
        + "/"
        + "coding.interval_list",
        fasta=config["resources"]["genome"],
    output:
        sbam=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["bams"]
            + "/"
            + "{patient}_split.out.bam"
        ),
    params:
        temporary_dir=config["TEMP_DIR"],
    threads: config["params"]["SplitNCigarReads"]["threads"]
    conda:
        "../envs/gatk.yml"
    resources:
        mem="32G",
        time="12:00:00",
        ncpus=4,
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["bam_cleaning"]
        + "/"
        + "{patient}.log",
    shell:
        """
        gatk --java-options "-Xmx30g -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" SplitNCigarReads -R {input.fasta} -I {input.bam} -O {output.sbam} \
        -L {input.intervals} --tmp-dir {params.temporary_dir}
        """
