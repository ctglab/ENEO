rule star_index:
    input:
        fasta=config["resources"]["genome"],
        gtf=config["resources"]["gtf"],
    output:
        directory(config["datadirs"]["index_folder"]),
    threads: config["params"]["STAR"]["threads"]
    params:
        uncompressed=lambda w, input: input[0].replace(".gz",""),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/star.yml"
    log:
        os.path.join(config["datadirs"]["logs"]["star_idx"], "star_idx.log"),
    resources:
        mem="60G",
        ncpus=8,
        runtime="360m",
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
        --genomeFastaFiles {params.uncompressed} --sjdbOverhang 100 --sjdbGTFfile {input.gtf}"""


rule salmon_gentrome:
    input:
        genome=config["resources"]["genome"],
        cdna=config["resources"]["transcriptome"],
    output:
        gentrome=temp(
            os.path.join(config["datadirs"]["salmon_idx"], "gentrome.fa.gz")
        ),
    log:
        os.path.join(config["datadirs"]["logs"]["salmon_quant"], "gentrome.log"),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/salmon.yml"
    shell:
        """
        cat {input.genome} {input.cdna} > {output.gentrome}
        """


rule salmon_idx:
    input:
        gentrome=os.path.join(config["datadirs"]["salmon_idx"], "gentrome.fa.gz"),
    output:
        out=os.path.join(config["datadirs"]["salmon_idx"], "ctable.bin"),
    threads: config["params"]["salmon"]["threads"]
    params:
        outdir=lambda w, output: os.path.dirname(os.path.abspath(output.out)),
        extra=config["params"]["salmon"]["extra"]["index"],
    resources:
        mem="40G",
        ncpus=8,
        runtime="240m",
    container:
        "docker://combinelab/salmon"
    conda:
        "../envs/salmon.yml"
    log:
        os.path.join(config["datadirs"]["logs"]["salmon_quant"], "index.log"),
    shell:
        """
        salmon index -t {input.gentrome} -i {params.outdir} {params.extra}
        """
