import os


rule vcfanno:
    input:
        toml_file=config["OUTPUT_FOLDER"] + "vcfanno.toml",
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "results",
            "variants",
            "variants.vcf.gz"
        ),
    output:
        vcf=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["VCF_out"],
                "{patient}_annot.vcf.gz"
            )
        ),
        index=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["VCF_out"],
                "{patient}_annot.vcf.gz.tbi"
            )
        ),
    params:
        vcfanno_binary=config["params"]["vcfanno"]["vcfanno_binary"],
        extra="--permissive-overlap",
        lua=config["params"]["vcfanno"]["vcfanno_lua"],
    threads: config["params"]["vcfanno"]["threads"]
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_vcfanno.log"
        ),
    resources:
        runtime="60m",
        ncpus=4,
        mem="16G",
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/vep.yml"
    shell:
        """
        {params.vcfanno_binary} -lua {params.lua} -p {threads} {params.extra} {input.toml_file} {input.vcf} |\
        bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule filtercalls:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot.vcf.gz"
        ),
        giab_intervals=config["resources"]["giab_intervals"],
    output:
        vcf=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["VCF_out"],
                "{patient}_DP_filt.vcf.gz"
            )
        ),
        index=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["VCF_out"],
                "{patient}_DP_filt.vcf.gz.tbi"
            )
        ),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_bcftools.log"
        ),
    resources:
        runtime="20m",
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
        toml_file=os.path.join(
            config["OUTPUT_FOLDER"],
            "vcfanno.toml"
        ),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/cyvcf2.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "toml.log"
        ),
    resources:
        runtime="20m",
        ncpus=1,
        mem="1G",
    shell:
        """
        python3 {input.toml_script} -y {input.config_main} -t {input.toml_template} -o {output.toml_file}
        """


rule germProb:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_DP_filt.vcf.gz"
        ),
        index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_DP_filt.vcf.gz.tbi"
        ),
        script=germProb_script,
    output:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz"
        ),
    container: "docker://danilotat/eneo"
    conda:
        "../envs/cyvcf2.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_germprob.log"
        ),
    resources:
        runtime="20m",
        ncpus=1,
        mem="4G",
    shell:
        """
        python3 {input.script} {input.vcf} {output.vcf}
        """


rule indexgermProb:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz"
        ),
    output:
        index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz.tbi"
        ),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/vep.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_idxgermProb.log"
        ),
    resources:
        runtime="20m",
        ncpus=1,
        mem="4G",
    shell:
        """
        tabix -p vcf {input.vcf}
        """
