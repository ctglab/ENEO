import os

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
        "docker://ctglabcnr/eneo"
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

rule MergeCalls:
    input:
        vcf1=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "deepvariant",
            "{patient}_deepvariant_FILT.vcf.gz",
        ),
        vcf2=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "strelka",
            "{patient}_strelka2_FILT.vcf.gz",
        ),
        ref_fasta=os.path.abspath(
            config["resources"]["genome"]
        ),
    output:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_merged_calls.vcf.gz",
        ),
        vcf_index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_merged_calls.vcf.gz.tbi",
        ),
    params:
        tmp_dv = os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_deepvariant_tmp.vcf.gz",
        ),
        tmp_strelka = os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],  
            "{patient}_strelka2_tmp.vcf.gz",
        ),
    container:
        "docker://ctglabcnr/eneo"
    conda:
        "../envs/samtools.yml"
    threads: config["params"]["samtools"]["threads"]
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["snv_calling"],
            "{patient}_merge_calls.log",
        ),
    resources:
        runtime="20m",
        ncpus=1,
        mem="8G",
    shell:
        """
        bcftools norm -f {input.ref_fasta} -m -both {input.vcf1} -Oz -o {params.tmp_dv}
        bcftools norm -f {input.ref_fasta} -m -both {input.vcf2} -Oz -o {params.tmp_strelka}
        tabix -p vcf {params.tmp_dv}
        tabix -p vcf {params.tmp_strelka}
        bcftools isec -n=2 -w2 {params.tmp_dv} {params.tmp_strelka} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        rm {params.tmp_dv} {params.tmp_dv}.tbi {params.tmp_strelka} {params.tmp_strelka}.tbi
        """


rule vcfanno:
    input:
        toml_file=os.path.join(
            config["OUTPUT_FOLDER"],
            "vcfanno.toml"),
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_merged_calls.vcf.gz",
        ),
        vcf_index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_merged_calls.vcf.gz.tbi",
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
        ncpus=1,
        mem="16G",
    container: "docker://ctglabcnr/eneo"
    conda: "../envs/vep.yml"
    shell:
        """
        {params.vcfanno_binary} -lua {params.lua} -p {threads} {params.extra} \
            {input.toml_file} {input.vcf} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule germProb:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot.vcf.gz"
        ),
        script=germProb_script,
    output:
        vcf=temp(
            os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz"
        )),
        vcf_idx=temp(
            os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz.tbi"
        )),
    container: "docker://ctglabcnr/eneo"
    conda: "../envs/cyvcf2.yml"
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
        tabix -p vcf {output.vcf}
        """
