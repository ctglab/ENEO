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
from concurrent.futures import ThreadPoolExecutor, as_completed

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
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel downloads/conversions.")
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
    conf["resources"][res_name] = os.path.abspath(path)
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

def process_chromosome_chunk(chrom_data, vcf_url, mapping_file, outfolder):
    """
    Worker function to process a single chromosome.
    """
    refseq_id, gencode_id = chrom_data
    chunk_filename = f"chunk_{gencode_id}.vcf.gz"
    chunk_path = os.path.join(outfolder, chunk_filename)
    
    if os.path.exists(chunk_path):
        logging.info(f"Chunk {gencode_id} already exists. Skipping.")
        return chunk_path, gencode_id
    logging.info(f"Processing chunk: {gencode_id} (Source: {refseq_id})")
    cmd = (
        f"bcftools view {vcf_url} -r {refseq_id} | "
        f"bcftools annotate --rename-chrs {mapping_file} | "
        f"bcftools sort -Oz -o {chunk_path}"
    )
    
    try:
        run_command(cmd, shell=True)
        # Index the chunk to ensure validity
        run_command(["tabix", "-p", "vcf", chunk_path])
        return chunk_path, gencode_id
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process chunk {gencode_id}: {e}")
        if os.path.exists(chunk_path):
            os.remove(chunk_path)
        raise e

def convert_notations(vcf_url, conv_from, conv_to, outfolder, max_workers=4):
    vcf_filename = vcf_url.split('/')[-1]
    out_vcf = os.path.join(outfolder, vcf_filename) # This will be the final Output
    
    if os.path.exists(out_vcf):
        logging.info(f"{out_vcf} already exists. Skipping conversion.")
        return out_vcf
        
    if "vcf" not in vcf_url:
        raise ValueError("Only VCF files supported.")
    if conv_from.lower() != "refseq" or conv_to.lower() != "gencode":
        raise NotImplementedError("Only refseq -> gencode is supported.")

    logging.info("Downloading chromosome alias table...")
    conv_df = pd.read_csv(
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz",
        sep="\t",
        compression="gzip",
        header=None,
        names=["source", "target", "source_type"]
    )
    mapping_df = conv_df[conv_df["source_type"] == "refseq"][["source", "target"]]
    mapping_file = os.path.join(outfolder, "refseq_dbsnp.tsv")
    mapping_df.to_csv(mapping_file, sep="\t", index=False, header=False)
    # just autosomal and sexual
    standard_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    target_map = mapping_df[mapping_df["target"].isin(standard_chroms)]
    # List of tuples: (NC_000001.11, chr1)
    chrom_tasks = list(zip(target_map["source"], target_map["target"]))
    logging.info(f"Starting parallel processing for {len(chrom_tasks)} chromosomes with {max_workers} workers.")
    chunk_files = {} # store paths for later
    # Run in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_chromosome_chunk, task, vcf_url, mapping_file, outfolder): task[1] 
            for task in chrom_tasks
        }
        for future in as_completed(futures):
            chrom_name = futures[future]
            try:
                path, _ = future.result()
                chunk_files[chrom_name] = path
            except Exception as e:
                logging.error(f"Exception for {chrom_name}: {e}")
                sys.exit(1)
    # merge chunks
    logging.info("Merging chunks into final VCF...")
    # Sort files based on standard chromosome order
    sorted_chunks = []
    for chrom in standard_chroms:
        if chrom in chunk_files:
            sorted_chunks.append(chunk_files[chrom])
    if not sorted_chunks:
        logging.error("No chunks were generated.")
        sys.exit(1)
    # Create a list file for concat
    list_file_path = os.path.join(outfolder, "chunks_list.txt")
    with open(list_file_path, "w") as f:
        for chunk in sorted_chunks:
            f.write(f"{chunk}\n")
    # concat chunks
    run_command(f"bcftools concat -f {list_file_path} -Oz -o {out_vcf}", shell=True)
    run_command(["tabix", "-p", "vcf", out_vcf])
    # drop chunks
    logging.info("Cleaning up temporary chunks...")
    for chunk in sorted_chunks:
        if os.path.exists(chunk):
            os.remove(chunk)
            if os.path.exists(chunk + ".tbi"):
                os.remove(chunk + ".tbi")
    if os.path.exists(list_file_path):
        os.remove(list_file_path)
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
    dict_file = f"{''.join(fasta_file.split('.')[:-1])}.dict" 
    fasta_file.replace(".fa", ".dict").replace(".fasta", ".dict")
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

def download_deepvariant_model_files(urls: list, outfolder: str):
    destpath = os.path.join(outfolder, "deepvariant_rna_model") 
    if not os.path.exists(destpath):
        os.makedirs(destpath, exist_ok=True)
    for url in urls:
        filename = url.split("/")[-1]
        if not os.path.exists(os.path.join(destpath, filename)):
            run_command(["wget", "-c", url, "-P", destpath])
    return destpath


def convert_REDI(bed_url, bed_output, drop_intermediate=True):
    if os.path.isfile(bed_output):
        logging.info(f"{bed_output} already exists.")
        return bed_output
    logging.info(f"Reading REDI portal file from {bed_url} and converting the annotation")
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
    
    # Ensure outfolder exists
    if not os.path.exists(outfolder):
        os.makedirs(outfolder, exist_ok=True)

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
                create_sequence_dictionary(path)
        elif ftype == "vcf":
            if name == "dbsnps":
                # conversion on-the-fly with parallelization
                path = convert_notations(res_entry['url'], "refseq", "gencode", outfolder, max_workers=args.workers)
                path = generate_allele_frequency(path, os.path.join(outfolder, f"{name}_withAF.vcf.gz"))
            else:
                path = download_resource(res_entry, outfolder, args.dry_run)
        elif ftype == "table":
            path = convert_REDI(res_entry["url"], os.path.join(outfolder, f"{name}.BED.gz"))
        elif ftype == "archive": 
            path = decompress_file(download_resource(res_entry, outfolder, args.dry_run))
        elif ftype == "model":
            path = download_deepvariant_model_files(res_entry['url'], outfolder)
        else:
            logging.warning(f"Unknown filetype for {name} as its {ftype}. Skipping.")
            continue
        conf = update_yaml(conf, name, path)
    with open(args.config, "w") as f:
        yaml.dump(conf, f)
    logging.info("Updated config written.")


if __name__ == "__main__":
    main(parse_arguments())