import os


rule vcfanno:
    input:
        toml_file=config["OUTPUT_FOLDER"] + "vcfanno.toml",
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow"
        + "/"
        + "results"
        + "/"
        + "variants"
        + "/"
        + "variants.vcf.gz",
    output:
        vcf=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_annot.vcf.gz"
        ),
        index=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_annot.vcf.gz.tbi"
        ),
    params:
        vcfanno_binary=config["params"]["vcfanno"]["vcfanno_binary"],
        extra="--permissive-overlap",
        lua=config["params"]["vcfanno"]["vcfanno_lua"],
    threads: config["params"]["vcfanno"]["threads"]
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["annotate_variants"]
        + "/"
        + "{patient}_vcfanno.log",
    resources:
        time="1:00:00",
        ncpus=4,
        mem="16G",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        {params.vcfanno_binary} -lua {params.lua} -p {threads} {params.extra} {input.toml_file} {input.vcf} |\
        bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule filtercalls:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot.vcf.gz",
        giab_intervals=config["resources"]["giab_intervals"],
    output:
        vcf=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_DP_filt.vcf.gz"
        ),
        index=temp(
            config["OUTPUT_FOLDER"]
            + config["datadirs"]["VCF_out"]
            + "/"
            + "{patient}_DP_filt.vcf.gz.tbi"
        ),
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["annotate_variants"]
        + "/"
        + "{patient}_bcftools.log",
    resources:
        time="0:20:00",
        ncpus=2,
        mem="8G",
    shell:
        """
        bcftools view -e "GT='mis'" {input.vcf} |\
         bcftools view -i "FILTER='PASS' & (DP > 5) & (FORMAT/AD[0:1] > 2)" --threads {threads} | bedtools intersect -header -v -a stdin -b {input.giab_intervals} \
         -sorted | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule createTOML:
    input:
        config_main=configpath,
        toml_template=config["params"]["vcfanno"]["vcfanno_toml"],
        toml_script=config["params"]["vcfanno"]["toml_script"],
    output:
        toml_file=config["OUTPUT_FOLDER"] + "vcfanno.toml",
    conda:
        "../envs/cyvcf2.yml"
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["annotate_variants"]
        + "/"
        + "toml.log",
    resources:
        time="0:20:00",
        ncpus=1,
        mem="1G",
    shell:
        """
        python3 {input.toml_script} -y {input.config_main} -t {input.toml_template} -o {output.toml_file}
        """


rule germProb:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_DP_filt.vcf.gz",
        index=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_DP_filt.vcf.gz.tbi",
        script=germProb_script,
    output:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz",
    conda:
        "../envs/cyvcf2.yml"
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["annotate_variants"]
        + "/"
        + "{patient}_germprob.log",
    resources:
        time="0:20:00",
        ncpus=1,
        mem="4G",
    shell:
        """
        python3 {input.script} {input.vcf} {output.vcf}
        """


rule indexgermProb:
    input:
        vcf=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz",
    output:
        index=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_annot_germProb.vcf.gz.tbi",
    conda:
        "../envs/samtools.yml"
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["annotate_variants"]
        + "/"
        + "{patient}_idxgermProb.log",
    resources:
        time="0:20:00",
        ncpus=1,
        mem="4G",
    shell:
        """
        tabix -p vcf {input.vcf}
        """
