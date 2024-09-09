"""
This script is used to read resources json file and the config file of the 
workflow, to download files and update accordingly the configuration file using
the updated resources.
"""

import yaml
import pandas as pd
import os
import subprocess
import sys
import json
import subprocess

cwd = os.path.abspath(__file__)

class ChromosomeConverter(object):
    """
    Basic class to control for chromosome notation and performing chromosome conversion
    based on the Ensembl notation.
    """
    def __init__(self) -> None:
        self.conv_table = self._get_table()
        self.possibilities = ["ensembl", "ucsc", "refseq"]

    def _get_table(self):
        conversion_table = pd.DataFrame()
        original_conv = pd.read_csv(
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz",
            compression="gzip", sep="\t", header=None, names=["from", "to", "how"])
        refseq_ucsc = original_conv.loc[original_conv["how"] == "refseq"]
        ensembl_ucsc = original_conv.loc[original_conv["how"].str.contains("ensembl")]
        refseq_ucsc["ensembl"] = refseq_ucsc["to"].map(ensembl_ucsc.set_index("to")["from"])
        refseq_ucsc["how"] = "refseq_to_ensembl"
        conversion_table = pd.concat([conversion_table, refseq_ucsc.loc[:, ["to", "ensembl", "how"]]], ignore_index=True).dropna(subset=["ensembl"])
        # now the same for the standard 
        ensembl_ucsc.rename(columns={"from": "ensembl"}, inplace=True)
        ensembl_ucsc["how"] = "ucsc_to_ensembl"
        conversion_table = pd.concat([conversion_table, ensembl_ucsc.loc[:,["to", "ensembl", "how"]]], ignore_index=True)
        conversion_table.rename(columns={"to": "from"}, inplace=True)
        return conversion_table
    
    @staticmethod
    def infer_notation_vcf(vcf_url: str):
        # this method infers the chromosome notation useful to return a from-to mask
        # it uses the ability of bcftools to read VCF on-the-fly
        cmd0 = f'bcftools view {vcf_url} | grep -v "#" | head -n 1 | cut -f1' 
        chromosome_notation = subprocess.run([cmd0], shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.strip()
        if chromosome_notation.startswith("chr"):
            return "ucsc"
        elif chromosome_notation.startswith("NC"):
            return "refseq"
        elif chromosome_notation.isnumeric():
            # even though this is prone to errors, this has to work because the VCF is sorted
            # so the first record is always the first chromosome
            return "ensembl"
        else:
            raise ValueError(f"{chromosome_notation} is not a recognized value.")


    def generate_txt(self, chr_from: str, chr_to: str, outfile: str):
        assert [x in self.possibilities for x in [chr_from, chr_to]]
        if any((x not in ["ensembl", "ucsc", "refseq"] for x in [chr_from, chr_to])):
            raise ValueError(f"Unknown value for {chr_from}")
        else:
            if not os.path.isfile(outfile):
                if not f"{chr_from}_to_{chr_to}" in self.conv_table["how"].unique():
                    raise ValueError(f"{chr_from}_to_{chr_to} is not in the conversion possibilities" )
                else:
                    subset = self.conv_table.loc[self.conv_table["how"] == f"{chr_from}_to_{chr_to}"]
                    subset.loc[:, ["from", "ensembl"]].to_csv(outfile, sep='\t', index=False, header=False)
                    print(f"{outfile} correctly created.")
            else:
                print(f"{outfile} already exists.")
    
    def convert_REDI(self, bed_input: str, bed_output: str, chr_to="ensembl"):
        chr_from = "ucsc"
        self.generate_txt(
            chr_from, chr_to, f"{chr_from}_to_{chr_to}.csv")
        conv_table = pd.read_csv(f"{chr_from}_to_{chr_to}.csv", sep='\t', header=None, names=["Region", "ensembl"])
        redi_df = pd.read_csv(bed_input, compression="gzip", sep='\t', usecols=[i for i in range(0,9)])
        # make it as a true BED, so converting with both start-end
        redi_df["Start"] = redi_df["Position"] - 1
        redi_df["End"] = redi_df["Start"] + 1
        redi_df = redi_df.loc[:, ["Region", "Start", "End"] + [x for x in redi_df.columns if x not in ["Region", "Position", "Start", "End"]]] 
        redi_df = redi_df.merge(conv_table, on="Region", how="left")
        print(f"{redi_df['ensembl'].isna().sum()} entries not converted")
        redi_df = redi_df.dropna(subset=["ensembl"]).drop(
            columns=["Region"]).rename(columns={"ensembl": "Region"})
        redi_df.loc[:, ["Region"] + [x for x in redi_df.columns if x != "Region"]].\
            to_csv(bed_output, index=False, header=False, sep='\t')
        # sort & compress + index
        cmd0 = f"sort -k1,1 -k2,2n -k3,3n {bed_output} | bgzip -c > {bed_output+'.gz'}"
        subprocess.run([cmd0], shell=True)
        subprocess.run(['tabix', '-p', 'bed', bed_output+".gz"])


class ResourceEntry(object):
    """
    This object handles the reading of the config file, check if file is present or 
    has to be downloaded. If a conversion has to be done, it will be performed. 
    Then it updates the relative configuration file accordingly.
    """
    def __init__(self, conf_file: dict, resources_entry, res_name: str, outfolder: str):
        self.conf_file = conf_file
        self.resources_entry = resources_entry
        self.res_name = res_name
        self.main_filename = self._get_main_filename()
        self.outfolder = outfolder
        self.res_type = self._get_restype()
        assert res_name in conf_file["resources"].keys()

    def _get_main_filename(self):
        if isinstance(self.resources_entry, dict):
            for file_url in self.resources_entry.values():
                if file_url.endswith(".gz"):
                    return file_url.split('/')[-1]
        else:
            return self.resources_entry.split('/')[-1] 
                

    def _get_restype(self):
        if isinstance(self.resources_entry, str):
            return "other"
        else:
            if "BED" in self.resources_entry.keys():
                return "BED"
            elif "VCF" in self.resources_entry.keys():
                return "VCF"
            else:
                raise ValueError("Unrecognized resource entry")

    def _download_stuff(self):
        if not os.path.isfile(self.conf_file["resources"][self.res_name]):
            if self.res_type in ["BED", "VCF"]:
                for record in self.resources_entry.values():
                    subprocess.run(['wget', '-c', record, '-P', self.outfolder])
            else:
                subprocess.run(['wget', '-c', self.resources_entry, '-P', self.outfolder])    
        else:
            print("File already downloaded")

def update_yaml(conf_main: str, resources: str, outfolder: str, chrom_converter: ChromosomeConverter):
    """
    This function does:
        - Reads configuration and resources file
        - Download stuff if needed
        - Update YAML accordingly 
    """
    with open(conf_main, "r") as conf_main_yaml:
        conf_main_yaml = yaml.load(conf_main_yaml, Loader=yaml.FullLoader)
    resources = json.load(open(resources, "r"))
    print(resources)
    for entry in conf_main_yaml["resources"]:
        if entry in resources.keys():
            resource_entry = ResourceEntry(conf_main_yaml, resources[entry], entry, outfolder)
            # for the VCF, read file on the fly.
            if resource_entry.res_type == "VCF":
                chr_notation = ChromosomeConverter.infer_notation_vcf(entry)
                if chr_notation == "ensembl":
                    # no conversion, just download
                    resource_entry._download_stuff()
                else:
                    # read and convert
                    cmd1 = f"bcftools view {entry} | bcftools annotate --rename-chrs {os.path.join(
                        outfolder, chr_notation)}_to_ensembl.csv | bgzip -c > {os.path.join(
                            outfolder, resource_entry.main_filename)}"
                    subprocess.run([cmd1], shell=True)
                    subprocess.run(["tabix", "-p", "vcf", os.path.join(outfolder, resource_entry.main_filename)])
                    # now update the yaml file
                    conf_main["resources"][entry] = os.path.join(outfolder, resource_entry.main_filename)

            # check chromosome notation. 
            # avoid indexes
            if resource_entry.res_type == "VCF":
                if resource_entry.main_filename.endswith(".gz"):
                    chr_notation = ChromosomeConverter.infer_notation(os.path.join(
                        outfolder, resource_entry.main_filename))
                    if chr_notation != "ensembl":
                        # we need to convert the file 
                        chrom_converter.convert_vcf(
                            vcf_input=os.path.join(outfolder, resource_entry.main_filename),
                            vcf_output=os.path.join(outfolder, resource_entry.main_filename.replace(
                            ".vcf.gz", "_converted.vcf.gz")))
                        conf_main_yaml["resources"][entry] = os.path.join(outfolder, resource_entry.main_filename.replace(".vcf.gz", "_converted.vcf.gz"))
                
            elif resource_entry.res_type == "BED":
                # we've got only the REDI portal file, which is not a properly BED file
                chrom_converter.convert_REDI(
                    bed_input=os.path.join(outfolder, resource_entry.main_filename),
                    bed_output=os.path.join(outfolder, resource_entry.main_filename.replace(
                        ".txt.gz", "_converted.BED")))
                conf_main_yaml["resources"][entry] = os.path.join(outfolder, resource_entry.main_filename.replace("txt.gz", "_converted.BED.gz"))
            else: 
                conf_main_yaml["resources"][entry] = os.path.join(outfolder, resource_entry.main_filename)
    # that's good, now we could write out the YAML 
    with open(conf_main, "w") as conf_main:
        yaml.dump(conf_main_yaml, conf_main)




if __name__ == "__main__":
    chr_converter = ChromosomeConverter()
    print(chr_converter.conv_table)
    # generate table
    chr_converter.generate_txt("ucsc", "ensembl", "ucsc_to_ensembl.csv")
    update_yaml(conf_main=sys.argv[1], resources=sys.argv[2], outfolder=os.getcwd(), chrom_converter=chr_converter)

