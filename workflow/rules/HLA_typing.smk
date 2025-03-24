# experimental execution of different paths based on .fastq.gz input or .bam

if execution_mode == "full":
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
        threads: config["params"]["t1k"]["threads"]
        resources:
            time="4:00:00",
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
elif execution_mode == "reduced":
    rule get_coords:
        input:
            genome=config["resources"]["genome"],
            gtf=config["resources"]["gtf"],
            idx=config["params"]["t1k"]["dat_file"],
        output:
            coords_dna=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "_dna_coord.fa"
            ),
            coords_rna=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "_rna_coord.fa"
            ),
        params:
            outdir=lambda w, output: os.path.dirname(os.path.abspath(output.coords_rna)),
        container:
            "docker://danilotat/eneo",
        resources:
            time="2:00:00",
            ncpus=4,
            mem="10G",
        log:
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["logs"]["t1k"],
                "coords.log"
            ),
        shell:
            """
            t1k-build.pl -f {input.genome} \
            -o {params.outdir} -d {input.idx} \
            -g {input.gtf}
            """
        
    rule genotype:
        input:
            unpack(get_bam),
            idx=config["resources"]["t1k_file"],
            coords=os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["HLA_typing"],
                "_rna_coord.fa"
            ),
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
        conda:
            "../envs/t1k.yml"
        threads: config["params"]["t1k"]["threads"]
        resources:
            time="4:00:00",
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
            run-t1k -b {input.bam} --preset hla \
            -f {input.idx} -c {input.coords} \
            -t {threads} -o {params.prefix} \
            --od {params.outdir} 
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
        log:
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["logs"]["t1k"],
                "{patient}_hla.log"
            ),
        resources:
            time="0:20:00",
            ncpus=2,
            mem="8G",
        shell:
            "python3 {input.hla_script} {input.genotype} > {output}"
