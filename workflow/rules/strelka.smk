rule Strelka_prep:
    input:
        cram=config["OUTPUT_FOLDER"]
        + config["datadirs"]["BQSR"]
        + "/"
        + "{patient}_recal.cram",
    output:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow"
        + "/"
        + "runWorkflow.py",
    params:
        ref_fasta=config["resources"]["genome"],
        regions=config["resources"]["intervals_coding"],
        runDir=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow",
    container:
        "docker://danilotat/strelka2"
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["snv_calling"]
        + "/"
        + "{patient}_strelka_prep.log",
    resources:
        time="0:20:00",
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
        script=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow"
        + "/"
        + "runWorkflow.py",
    output:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow"
        + "/results/variants/"
        + "variants.vcf.gz",
        txt=config["OUTPUT_FOLDER"]
        + config["datadirs"]["VCF_out"]
        + "/"
        + "{patient}_workflow"
        + "/results/"
        + "checkpoint.txt",
    params:
        threads=config["params"]["strelka2"]["threads"],
    container:
        "docker://danilotat/strelka2"
    log:
        config["OUTPUT_FOLDER"]
        + config["datadirs"]["logs"]["snv_calling"]
        + "/"
        + "{patient}_calling.log",
    resources:
        time="4:00:00",
        ncpus=2,
        mem="16G",
    shell:
        """
        {input.script} -m local -j {params.threads}
        echo "finished" > {output.txt}
        """
