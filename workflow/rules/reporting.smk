import os 
from pathlib import Path

rule multiqc:
    input:
        unpack(get_multiqc_inputs),
    output:
        html=os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["qc_reports"],
            "multiqc_report.html"
        ),
        data=directory(
            os.path.join(
                config["OUTPUT_FOLDER"],
                config["datadirs"]["qc_reports"],
                "multiqc_data"
            )
        ),
    params:
        outdir=lambda wc, output: Path(
            output.data).parent.absolute(),
        title="ENEO QC Report",
    container:
        "docker://ewels/multiqc:latest"
    conda:
        "../envs/multiqc.yml"
    resources:
        mem="8G",
        runtime="60m",
        ncpus=1,
    log:
        os.path.join(
            config["OUTPUT_FOLDER"],
            config["datadirs"]["logs"]["trimming"],
            "multiqc.log"
        ),
    shell:
        """
        multiqc {input} \
            --outdir {params.outdir} \
            --title "{params.title}" \
            --force \
            2>&1 | tee {log}
        """
