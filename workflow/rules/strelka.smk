import os

rule Strelka_prep:
    input:
        bam=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.bam"
        ),
    output:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "runWorkflow.py",
        ),
    params:
        ref_fasta=os.path.abspath(config["resources"]["genome"]),
        regions=config["resources"]["intervals_coding"],
        runDir=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
        ),
    conda:
        "../envs/strelka2.yml"
    container: "docker://ctglabcnr/strelka2"
    
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
        --bam {input.bam} \
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
        bam=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.bam"
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
    container: "docker://swantonlab/strelka2"
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

rule SelectStrelka2Calls:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_workflow",
            "results/variants/variants.vcf.gz",
        ),
        giab_intervals=os.path.abspath(
            config["resources"]["giab_intervals"]
        ),
    output:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_strelka2_FILT.vcf.gz",
        ),
        vcf_index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "{patient}_strelka2_FILT.vcf.gz.tbi",
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
            "{patient}_select_strelka2_calls.log",
        ),
    resources:
        runtime="20m",
        ncpus=2,
        mem="8G",
    shell:
        """
        bcftools view -e "GT='mis'" {input.vcf} |\
         bcftools view -i "FILTER='PASS' & (DP > 5) & (FORMAT/AD[0:1] > 2)" --threads {threads} | \
         bedtools intersect -header -v -a stdin -b {input.giab_intervals} -sorted | \
         bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """