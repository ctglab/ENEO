rule star_index:
    input:
        fasta=config["resources"]["genome"],
        gtf=config["resources"]["gtf"],
    output:
        directory(config["datadirs"]["index_folder"]),
    threads: config["params"]["STAR"]["threads"]
    container:
        "docker://danilotat/eneo"
    log:
        config["datadirs"]["logs"]["star_idx"]
        + "/"
        + "star_idx.log",
    resources:
        mem="60G",
        ncpus=8,
        time="6:00:00",
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
        --genomeFastaFiles genome --sjdbOverhang 100 --sjdbGTFfile {input.gtf}"""


rule salmon_gentrome:
    input:
        genome=config["resources"]["genome"],
        cdna=config["resources"]["transcriptome"],
    output:
        temp(
            config["datadirs"]["salmon_idx"] 
            + "/" 
            + "gentrome.fa.gz"),
    log:
        config["datadirs"]["logs"]["salmon_quant"] 
        + "/" 
        + "gentrome.log",
    conda:
        "../envs/salmon_new.yml"
    shell:
        "cat {input.genome} {input.cdna} | gzip > {output}"


rule salmon_idx:
    input:
        gentrome=config["datadirs"]["salmon_idx"] 
        + "/" 
        + "gentrome.fa.gz",
    output:
        out=config["datadirs"]["salmon_idx"] 
        + "/" 
        + "ctable.bin",
    threads: config["params"]["salmon"]["threads"]
    params:
        outdir=lambda w, output: os.path.dirname(os.path.abspath(output.out)),
        extra=config["params"]["salmon"]["extra"]["index"],
    resources:
        mem="40G",
        ncpus=8,
        time="4:00:00",
    conda:
        "../envs/salmon_new.yml"
    log:
        config["datadirs"]["logs"]["salmon_quant"] 
        + "/"
        + "index.log",
    shell:
        """
        salmon index -t {input.gentrome} -i {params.outdir} \
        --decoys {input.decoys} {params.extra}
        """
