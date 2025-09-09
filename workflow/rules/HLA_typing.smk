# experimental execution of different paths based on .fastq.gz input or .bam

rule genotype:
    input:
        unpack(get_fastq),
        idx=config["resources"]["t1k_file"],
    output:
        temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "{patient}_aligned_1.fa"
            )
        ),
        temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "{patient}_aligned_2.fa"
            )
        ),
        temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "{patient}_allele.tsv"
            )
        ),
        temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "{patient}_allele.vcf"
            )
        ),
        temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "{patient}_candidate_1.fq"
            )
        ),
        temp(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "{patient}_candidate_2.fq"
            )
        ),
        hla=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["HLA_typing"],
            "{patient}_genotype.tsv"
        ),
    params:
        prefix="{patient}",
        outdir=lambda w, output: os.path.dirname(os.path.abspath(output.hla)),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/t1k.yml"
    threads: config["params"]["t1k"]["threads"]
    resources:
        runtime="240m",
        ncpus=4,
        mem="32G",
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["t1k"],
            "{patient}.log"
        ),
    shell:
        """
        run-t1k -1 {input.r1} -2 {input.r2} --preset hla \
        -f {input.idx} -t {threads} -o {params.prefix} --od {params.outdir}
        """
        
rule extract_hla:
    input:
        genotype=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["HLA_typing"],
            "{patient}_genotype.tsv"
        ),
        hla_script=config["resources"]["hla_script"],
    output:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["HLA_typing"],
            "{patient}_allele_input_pvacseq.csv"
        ),
    container:
        "docker://danilotat/eneo"
    conda:
        "../envs/cyvcf2.yml"
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["t1k"],
            "{patient}_hla.log"
        ),
    resources:
        runtime="20m",
        ncpus=2,
        mem="8G",
    shell:
        "python3 {input.hla_script} {input.genotype} > {output}"
