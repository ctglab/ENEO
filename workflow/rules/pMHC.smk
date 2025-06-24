import os

rule pMHCpeptides:
    input:
        vcf=os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["VCF_out"], "{patient}.vep.vcf.gz"),
        vcf_idx=os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["VCF_out"], "{patient}.vep.vcf.gz.tbi"),
        hla=os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["HLA_typing"], "{patient}_allele_input_pvacseq.csv"),
        launcher=config["params"]["pMHC"]["netmhcpan_launcher_script"],
    params:
        patname="{patient}",
        outfolder=lambda w, output: os.path.dirname(os.path.abspath(output.out)),
        threads=config["params"]["pMHC"]["threads"],
    output:
        out=os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["peptides"], "{patient}", "{patient}.epitopes.csv"),
    container:
        "docker://danilotat/eneo"
    log:
        os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["logs"]["pMHC"], "{patient}.log"),
    resources:
        runtime="120m",
        ncpus=4,
        mem="8G",
        tmpdir=config["TEMP_DIR"],
    shell:
        """
        python3 {input.launcher} --vcf {input.vcf} -a {input.hla} -p {params.patname} -o {params.outfolder} -n {params.threads}
        """
    
rule filter_peptides:
    input:
        peptides=os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["peptides"], "{patient}", "{patient}.epitopes.csv"),
        calibration_frame=config["params"]["pMHC"]["calibration_frame"],
        hla_ligand_atlas=config["params"]["pMHC"]["hla_ligand_atlas"],
        script=config["params"]["pMHC"]["filter_peptides_script"],
    params:
        min_length=config["params"]["pMHC"]["min_length"],
        max_length=config["params"]["pMHC"]["max_length"],
        germProb=config["params"]["pMHC"]["germProb"],
    output:
        out=os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["peptides"], "{patient}.filtered_epitopes.csv"),
    log:
        os.path.join(config["OUTPUT_FOLDER"], config["datadirs"]["logs"]["pMHC"], "{patient}_filt.log"),
    resources:
        runtime="60m",
        ncpus=2,
        mem="2G",
    container:
        "docker://danilotat/eneo"
    shell:
        """
        python3 {input.script} -i {input.peptides} -c {input.calibration_frame} -a {input.hla_ligand_atlas} -o {output.out} -min {params.min_length} -max {params.max_length} -g {params.germProb}
        """
