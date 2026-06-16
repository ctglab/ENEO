rule star_index:
    input:
        fasta=os.path.abspath(
            config["resources"]["genome"]).rstrip(".gz"),
        gtf=config["resources"]["gtf"],
    output:
        directory(config["datadirs"]["index_folder"]),
    threads: config["params"]["STAR"]["threads"]
    params:
    container:
        "docker://ctglabcnr/star"
    conda:
        "../envs/star.yml"
    log:
        os.path.join(config["datadirs"]["logs"]["star_idx"], "star_idx.log"),
    resources:
        mem="60G",
        ncpus=1,
        runtime="360m",
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
        --genomeFastaFiles {input.fasta} --sjdbOverhang 100 --sjdbGTFfile {input.gtf}
        """


rule salmon_gentrome:
    input:
        genome=config["resources"]["genome"],
        cdna=config["resources"]["transcriptome"],
    output:
        gentrome=temp(
            os.path.join(config["datadirs"]["salmon_idx"], "gentrome.fa.gz")
        ),
        decoys=temp(
            os.path.join(config["datadirs"]["salmon_idx"], "decoys.txt")
        ),
    log:
        os.path.join(config["datadirs"]["logs"]["salmon_quant"], "gentrome.log"),
    container:
        "docker://ctglabcnr/eneo"
    conda:
        "../envs/salmon.yml"
    shell:
        """
        # Decoy names are the genome sequence ids (first token of each header).
        grep '^>' {input.genome} | sed 's/^>//' | cut -d ' ' -f1 > {output.decoys}
        # Decoy-aware gentrome: transcriptome targets first, genome (decoys) last.
        # GENCODE cdna headers are stripped to the bare transcript id so salmon's
        # target names match the tx2gene mapping used downstream.
        cat <(sed -E 's/^>([^|]+).*/>\\1/' {input.cdna}) {input.genome} > {output.gentrome}
        """


rule salmon_idx:
    input:
        gentrome=os.path.join(config["datadirs"]["salmon_idx"], "gentrome.fa.gz"),
        decoys=os.path.join(config["datadirs"]["salmon_idx"], "decoys.txt"),
    output:
        out=os.path.join(config["datadirs"]["salmon_idx"], "ctable.bin"),
    threads: config["params"]["salmon"]["threads"]
    params:
        outdir=lambda w, output: os.path.dirname(os.path.abspath(output.out)),
        extra=config["params"]["salmon"]["extra"]["index"],
    resources:
        mem="40G",
        ncpus=1,
        runtime="240m",
    container:
        "docker://combinelab/salmon"
    conda:
        "../envs/salmon.yml"
    log:
        os.path.join(config["datadirs"]["logs"]["salmon_quant"], "index.log"),
    shell:
        """
        salmon index -t {input.gentrome} -d {input.decoys} -i {params.outdir} -p {threads} {params.extra}
        """
