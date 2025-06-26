#!/usr/bin/env python3

import os
import sys
import yaml
import json
import argparse
import logging
import subprocess
from pathlib import Path
import pandas as pd
from rich.logging import RichHandler

# Setup logging
logging.basicConfig(
    level="DEBUG",
    format="%(message)s",
    datefmt="[%Y-%m-%d %H:%M:%S]",
    handlers=[RichHandler(markup=False, rich_tracebacks=True)]
)
cwd = Path(__file__).parent
config_folder = cwd.parent / "config"

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download and setup ENEO resources.")
    parser.add_argument("-c", "--config", type=str, default=config_folder / "config_main.yaml", help="Path to config.")
    parser.add_argument("-r", "--resources", type=str, default=cwd / "resources.json", help="Path to resources JSON.")
    parser.add_argument("-o", "--outfolder", type=str, required=True, help="Output folder.")
    parser.add_argument("--dry-run", action="store_true", help="Only simulate.")
    return parser.parse_args()


def read_yaml(path):
    with open(path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    for key in ["OUTPUT_FOLDER", "TEMP_DIR"]:
        if not data.get(key):
            logging.warning(f"{key} is not set in the configuration.")
    return data


def read_json(path):
    with open(path) as f:
        return json.load(f)

def update_yaml(conf, res_name, path):
    if res_name not in conf.get("resources", {}):
        logging.error(f"{res_name} not in config.")
        sys.exit(1)

    logging.info(f"Updating config: {res_name}")
    conf["resources"][res_name] = path
    return conf

def run_command(cmd, shell=False):
    logging.debug(f"Running command: {cmd}")
    subprocess.run(cmd, check=True, shell=shell)

def download_resource(entry, outfolder, dry=False):
    filename = entry["url"].split("/")[-1]
    dest = os.path.join(outfolder, filename)
    if not os.path.isfile(dest):
        logging.info(f"Downloading {filename}")
        if not dry:
            run_command(["wget", "-c", entry["url"], "-P", outfolder])
            if entry["url"].endswith('vcf.gz'):
                # download also the index
                index_url = entry["url"] + ".tbi"
                run_command(["wget", "-c", index_url, "-P", outfolder])
    else:
        logging.info(f"{filename} already exists.")
    return dest


def decompress_file(path):
    if path.endswith((".gtf.gz", ".fasta.gz", ".fa.gz")):
        target = path.replace(".gz", "")
        if not os.path.isfile(target):
            logging.info(f"Decompressing {path}")
            run_command(["gzip", "-d", path])
        return target
    if path.endswith(".tar.gz"):
        target_dir = path.replace(".tar.gz", "")
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir, exist_ok=True)
            logging.info(f"Extracting archive {path}")
            run_command(["tar", "-xzvf", path, "-C", target_dir])
        return target_dir
    return path


def convert_notations(vcf_url, conv_from, conv_to, outfolder):
    vcf_path = os.path.join(outfolder, vcf_url.split('/')[-1])
    if os.path.exists(os.path.join(outfolder, vcf_path)):
        logging.info(f"{vcf_path} already exists. Skipping conversion.")
        return vcf_path
    if "vcf" not in vcf_url:
        raise ValueError("Only VCF files supported.")
    if conv_from.lower() != "refseq" or conv_to.lower() != "gencode":
        raise NotImplementedError("Only refseq -> gencode is supported.")
    conv_df = pd.read_csv(
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz",
        sep="\t",
        compression="gzip",
        header=None,
        names=["source", "target", "source_type"]
    )
    mapping = conv_df[conv_df["source_type"] == "refseq"][["source", "target"]]
    mapping_file = os.path.join(outfolder, "refseq_dbsnp.tsv")
    mapping.to_csv(mapping_file, sep="\t", index=False, header=False)
    out_vcf = os.path.join(outfolder, os.path.basename(vcf_path))
    run_command(
        f"bcftools annotate --rename-chrs {mapping_file} {vcf_url} | bcftools sort -Oz -o {out_vcf}",
        shell=True)
    run_command(["tabix", "-p", "vcf", out_vcf])
    return out_vcf

def generate_allele_frequency(input_vcf, output_vcf):
    logging.info(f"Adding allele frequency to {input_vcf}")
    tags = "'INFO/AN:1=int(smpl_sum(AN))','INFO/AC:1=int(smpl_sum(AC))','INFO/AF:1=float(AC/AN)'"
    run_command(f"bcftools +fill-tags {input_vcf} -Oz -o {output_vcf} -- -t {tags}", shell=True)
    run_command(["tabix", "-p", "vcf", output_vcf])
    return output_vcf

def create_sequence_dictionary(fasta_file):
    """
    Create index and sequence dictionary for a FASTA file using samtools 
    """
    dict_file = fasta_file.replace(".fa", ".dict").replace(".fasta", ".dict")
    index_file = fasta_file + ".fai"
    for file in [dict_file, index_file]:
        if os.path.isfile(file):
            logging.info(f"{file} already exists. Skipping creation.")
        else:
            logging.info(f"Creating index for {fasta_file}")
            if file == index_file:
                run_command(["samtools", "faidx", fasta_file])
            else:
                logging.info(f"Creating sequence dictionary for {fasta_file}")
                run_command(["samtools", "dict", fasta_file, "-o", dict_file])

def convert_REDI(bed_url, bed_output, drop_intermediate=True):
    if os.path.isfile(bed_output):
        logging.info(f"{bed_output} already exists.")
        return bed_output

    redi_df = pd.read_csv(bed_url, compression="gzip", sep="\t", usecols=range(9))
    redi_df["Start"] = redi_df["Position"] - 1
    redi_df["End"] = redi_df["Start"] + 1
    bed_df = redi_df[["Region", "Start", "End"] + [col for col in redi_df.columns if col not in ["Region", "Position", "Start", "End"]]]
    bed_df.sort_values(["Region", "Start", "End"], inplace=True)
    tmp_output = bed_output.replace(".gz", "")
    bed_df.to_csv(tmp_output, sep="\t", index=False, header=False)
    run_command(f"cat {tmp_output} | bgzip -c > {bed_output}", shell=True)
    run_command(["tabix", "-p", "bed", bed_output])
    if drop_intermediate and os.path.isfile(tmp_output):
        os.remove(tmp_output)
    return bed_output


def main(args):
    conf = read_yaml(args.config)
    resources = read_json(args.resources)
    outfolder = args.outfolder
    for name, existing_path in conf.get("resources", {}).items():
        if name not in resources and not os.path.isfile(existing_path):
            logging.error(f"{name} missing in resources and not in repo.")
            continue
        if os.path.isfile(existing_path):
            logging.info(f"{name} already exists. Skipping.")
            continue
        res_entry = resources.get(name)
        if not res_entry:
            continue
        ftype = res_entry["filetype"].lower()
        if ftype in {"fasta", "gtf"}:
            path = decompress_file(download_resource(res_entry, outfolder, args.dry_run))
            if name == "genome":
                # create index and sequence dictionary for FASTA
                create_sequence_dictionary(path)
        elif ftype == "vcf":
            if name == "dbsnps":
                # conversion on-the-fly instead of downloading
                path = convert_notations(res_entry['url'], "refseq", "gencode", outfolder)
                path = generate_allele_frequency(path, os.path.join(outfolder, f"{name}_withAF.vcf.gz"))
            else:
                path = download_resource(res_entry, outfolder, args.dry_run)
        elif ftype == "table":
            path = convert_REDI(res_entry["url"], os.path.join(outfolder, f"{name}.BED.gz"))
        elif ftype == "archive":
            path = decompress_file(download_resource(res_entry, outfolder, args.dry_run))
        else:
            logging.warning(f"Unknown filetype for {name} as its {ftype}. Skipping.")
            continue
        conf = update_yaml(conf, name, path)
    with open(args.config, "w") as f:
        yaml.dump(conf, f)
    logging.info("Updated config written.")


if __name__ == "__main__":
    main(parse_arguments())
