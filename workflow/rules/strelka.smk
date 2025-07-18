import os

rule Strelka_prep:
    input:
        cram=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.cram",
        ),
    output:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "runWorkflow.py",
        ),
    params:
        ref_fasta=config["resources"]["genome"],
        regions=config["resources"]["intervals_coding"],
        runDir=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
        ),
    container:
        "docker://danilotat/strelka2",
    conda:
        "../envs/strelka2.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["snv_calling"],
            "{patient}_strelka_prep.log",
        ),
    resources:
        runtime="20m",
        ncpus=1,
        mem="4G",
    shell:
        """
        configureStrelkaGermlineWorkflow.py \
        --bam {input.cram} \
        --rna \
        --referenceFasta {params.ref_fasta} \
        --callRegions {params.regions} \
        --runDir {params.runDir} \
        --reportEVSFeatures
        """


rule Strelka2:
    input:
        script=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "runWorkflow.py",
        ),
    output:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "results/variants/variants.vcf.gz",
        ),
        txt=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "results/checkpoint.txt",
        ),
    params:
        threads=config["params"]["strelka2"]["threads"],
    container:
        "docker://danilotat/strelka2",
    conda:
        "../envs/strelka2.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["snv_calling"],
            "{patient}_calling.log",
        ),
    resources:
        runtime="240m",
        ncpus=2,
        mem="16G",
    shell:
        """
        {input.script} -m local -j {params.threads}
        echo "finished" > {output.txt}
        """
