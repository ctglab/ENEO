import os 


rule annotate_variants:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz"
        ),
        vcf_idx=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_annot_germProb.vcf.gz.tbi"
        ),
        plugin_wt=config["params"]["vep"]["extra"]["plugins"]["Wildtype"],
        plugin_fs=config["params"]["vep"]["extra"]["plugins"]["Frameshift"],
        cache=config["resources"]["vep_cache"],
    output:
        vcfout=temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["VCF_out"],
                "{patient}.annotated.vcf"
                )
        ),
    params:
        assembly=config["params"]["vep"]["extra"]["assembly"],
        filtering=config["params"]["vep"]["extra"]["filtering"],
        plugin_dir=lambda wc, input: os.path.dirname(input.plugin_wt),
    container:
        "docker://danilotat/eneo",
    conda:
        "../envs/vep.yml",
    resources:
        mem="6G",
        time="2:00:00",
        ncpus=2,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}.log"
        ),
    shell:
        """ 
        vep --input_file {input.vcf} \
        --output_file {output.vcfout} \
        --format vcf --vcf --symbol --terms SO --tsl \
        {params.filtering} \
        --plugin Frameshift --plugin Wildtype \
        --dir_plugins {params.plugin_dir} \
        --force_overwrite \
        --assembly {params.assembly} \
        --offline --cache --dir_cache {input.cache} 
        """

rule compress_annotated_vcf:
    input:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}.annotated.vcf"
        ),
    output:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}.vep.vcf.gz"
        ),
        vcf_idx=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}.vep.vcf.gz.tbi"
        ),
    container:
        "docker://danilotat/eneo",
    conda:
        "../envs/vep.yml",
    resources:
        mem="6G",
        time="1:00:00",
        ncpus=2,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_compress.log"
        ),
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """
        
