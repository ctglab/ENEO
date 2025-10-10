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
        # branching for do annotation online in CI mode. 
        cache=branch(
            config['execution_mode'] == 'CI',
            then="",
            otherwise=config["resources"]["vep_cache"],
        ),
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
        cache=lambda wc, input: "--database" if config['execution_mode'] == 'CI' else f"--offline --cache --dir_cache {os.path.dirname(input.cache)}",
        plugin_dir=lambda wc, input: os.path.dirname(input.plugin_wt),
    container:
        "docker://ensemblorg/ensembl-vep:release_105.0",
    conda:
        "../envs/vep.yml",
    resources:
        mem="6G",
        runtime="120m",
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
        {params.cache} \
        """

rule compress_annotated_vcf:
    input:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}.annotated.vcf"
        ),
    output:
        vcf=temp(
            os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}.vep.vcf.gz"
        )),
        vcf_idx=temp(
            os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}.vep.vcf.gz.tbi"
        )),
    container:
        "docker://ctglabcnr/eneo",
    conda:
        "../envs/vep.yml",
    resources:
        mem="6G",
        runtime="60m",
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
        
rule rna_errors:
    input:
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
    output:
        vcfout=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_final.vcf.gz"
        ),
        vcfout_idx=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_final.vcf.gz.tbi"
        ),
    params:
        script=config["resources"]["rna_errors_script"],
        reference=os.path.abspath(
            config["resources"]["genome"]).rstrip(".gz"),
        PoN=config["resources"]["PoN"],
        gtf=config["resources"]["gtf"],
        giab=config["resources"]["giab_intervals"],
        patID="{patient}",
    container:
        "docker://ctglabcnr/eneo",
    conda:
        "../envs/cyvcf2.yml",
    resources:
        mem="6G",
        runtime="60m",
        ncpus=2,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_rna_errors.log"
        ),
    shell:
        """
        python {params.script} \
            --vcf {input.vcf} \
            --ref {params.reference} \
            --giab {params.giab} \
            --gtf {params.gtf} \
            --pon {params.PoN} \
            --pat {params.patID} \
            --out {output.vcfout}
        tabix -p vcf {output.vcfout}
        """

rule passonly:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_final.vcf.gz"
        ),
        vcf_idx=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_final.vcf.gz.tbi"
        ),
    output:
        vcfout=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_final_passonly.vcf.gz"
        ),
        vcfout_idx=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_final_passonly.vcf.gz.tbi"
        ),
    container:
        "docker://ctglabcnr/eneo",
    conda:
        "../envs/cyvcf2.yml",
    resources:
        mem="6G",
        runtime="60m",
        ncpus=2,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["annotate_variants"],
            "{patient}_passonly.log"
        ),
    shell:
        """
        bcftools view -f .,PASS {input.vcf} -Oz -o {output.vcfout}
        tabix -p vcf {output.vcfout}
        """