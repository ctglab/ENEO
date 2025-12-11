import os

rule DeepVariant:
    input:
        bam=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["BQSR"],
            "{patient}_recal.bam"
        ),
        ref_fasta=os.path.abspath(config["resources"]["genome"]),
        regions=config["resources"]["intervals_coding"],
    output:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "deepvariant",
            "{patient}_deepvariant.vcf.gz",
        ),
        vcf_index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "deepvariant",
            "{patient}_deepvariant.vcf.gz.tbi",
        ),
    params:
        rna_model=os.path.join(
            os.path.abspath(config["resources"]["deepvariant_rna_model"]), "model.ckpt"
        ),
        uncompressed_regions=lambda wc, input:
            input.regions.replace('.gz',''),
        tmp_dir=os.path.join(
            config["OUTPUT_FOLDER"],
            "tmp",
            "{patient}_deepvariant_tmp",
        ),
        sample_name="{patient}",
        threads=config["params"]["deepvariant"]["threads"],
        extra=config["params"]["deepvariant"]["extra"],
    container:
        "docker://google/deepvariant:1.4.0",
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["snv_calling"],
            "{patient}_deepvariant.log",
        ),
    resources:
        runtime="480m",
        ncpus=config["params"]["deepvariant"]["threads"],
        mem="16G",
    shell:
        """
        mkdir -p {params.tmp_dir}
        run_deepvariant \
        --model_type=WES \
        --customized_model={params.rna_model} \
        --ref={input.ref_fasta} \
        --reads={input.bam} \
        --output_vcf={output.vcf} \
        --num_shards={params.threads} \
        --regions={params.uncompressed_regions} \
        --make_examples_extra_args={params.extra} \
        --intermediate_results_dir {params.tmp_dir} 
        rm -rf {params.tmp_dir}
        """

rule SelectDeepVariantCalls:
    input:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "deepvariant",
            "{patient}_deepvariant.vcf.gz",
        ),
        giab_intervals=os.path.abspath(
            config["resources"]["giab_intervals"]
        ),
    output:
        vcf=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "deepvariant",
            "{patient}_deepvariant_FILT.vcf.gz",
        ),
        vcf_index=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["VCF_out"],
            "deepvariant",
            "{patient}_deepvariant_FILT.vcf.gz.tbi",
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
            "{patient}_select_deepvariant_calls.log",
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
       
