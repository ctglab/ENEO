import os 

rule salmon_quantification:
    input:
        unpack(get_fastq),
        index=os.path.join(config["datadirs"]["salmon_idx"], "ctable.bin"),
    output:
        quant=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["salmon_quant"],
            "{patient}",
            "quant.sf",
        ),
    params:
        index=lambda wc, input: os.path.dirname(os.path.abspath(input.index)),
        libtype=config["params"]["salmon"]["extra"]["libtype"],
        zip_ext=config["params"]["salmon"]["extra"]["zip_ext"],
        extra=config["params"]["salmon"]["extra"]["extra"],
        outdir=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["salmon_quant"],
            "{patient}",
        ),
    threads: config["params"]["salmon"]["threads"]
    resources:
        time="1:00:00",
        ncpus=4,
        mem="32G",
    conda:
        "../envs/salmon_new.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["salmon_quant"],
            "{patient}.log",
        ),
    shell:
        """
        salmon quant -l {params.libtype} -i {params.index} -1 {input.r1} -2 {input.r2} -p {threads} -o {params.outdir}
        """


rule export_quantification:
    input:
        quant=expand(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["salmon_quant"],
                "{patient}",
                "quant.sf",
            ),
            patient=patients,
        ),
        idx=os.path.join(config["datadirs"]["salmon_idx"], "ctable.bin"),
        index=lambda wc, input: os.path.dirname(os.path.abspath(input.index)),
        cdna_fasta=config["resources"]["transcriptome"],
        annotation=config["resources"]["gtf"],
    output:
        transcript=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["expression"],
            "transcript_expression.tsv",
        ),
        gene=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["expression"],
            "gene_expression.tsv",
        ),
    params:
        outfolder=lambda w, output: os.path.dirname(os.path.abspath(output.gene)),
        patients=expand("{patient}", patient=patients),
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["export_quant"],
            "export_quantification.log",
        ),
    resources:
        time="0:30:00",
        ncpus=2,
        mem="8G",
    conda:
        "../envs/merge_salmon_quant.yml"
    script:
        "../scripts/merge_salmon_quantification.R"
