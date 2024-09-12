"""
This script is intended to be used as a one-shot resource configuration worker to setup all the needed files to run ENEO. As everything is formatted following the Ensembl annotation, whenever needed each file will be converted to match requirements. 
Report any issue on Github at:
    github.com/ctglab/ENEO/issues
"""

import yaml
import pandas as pd
import os
import subprocess
import sys
import json
import subprocess
import cyvcf2

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
    

    def convert_REDI(self, bed_url: str, bed_output: str, chr_to="ensembl", drop_intermediate=True):
        chr_from = "ucsc"
        self.generate_txt(
            chr_from, chr_to, f"{chr_from}_to_{chr_to}.csv")
        conv_table = pd.read_csv(f"{chr_from}_to_{chr_to}.csv", sep='\t', header=None, names=["Region", "ensembl"])
        redi_df = pd.read_csv(bed_url, compression="gzip", sep='\t', usecols=[i for i in range(0,9)])
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
        if drop_intermediate:
            os.remove(bed_output)

class ResourceEntry(object):
    """
    This object handles the reading of the config file, check if file is present or 
    has to be downloaded. If a conversion has to be done, it will be performed. 
    Then it updates the relative configuration file accordingly.
    """
    def __init__(self, conf_file: dict, resources_entry: dict, res_name: str, outfolder: str):
        self.conf_file = conf_file
        self.resources_entry = resources_entry
        self.res_name = res_name
        self.main_filename = self._get_main_filename()
        self.outfolder = outfolder
        self.downloaded = self._is_downloaded()
        self.res_type = self._get_restype()
        assert res_name in conf_file["resources"].keys()

    def _get_main_filename(self):
        return self.resources_entry["url"].split('/')[-1]
                
    def _get_restype(self):
        try:
            filetype = self.resources_entry["filetype"]
            return filetype
        except KeyError:
            print(f"Malformed entry for {self.res_name}. Check the resource file")

    def _is_downloaded(self):
        if os.path.isfile(os.path.join(self.outfolder, self.main_filename)):
            return True
        else:
            return False

    def _download_stuff(self):
        if not self.downloaded:
            print(f"Downloading {self.res_name}")
            subprocess.run(['wget', '-c', self.resources_entry["url"], '-P', self.outfolder])
            self.downloaded = True
        else:
            print(f"{self.res_name} already downloaded.")
            pass            


    def generate_allele_frequency(self):
        """
        This function edits the dbSNP VCF to generate three distinct INFO fields
        - AN: Alleles number
        - AC: Alleles count
        - AF: Alleles frequency
        """
        def calc_freq(variant: cyvcf2.Variant):
            ANs = variant.format(field='AN')
            ACs = variant.format(field='AC')
            # the AN and AC fields are lists of strings
            # I'll take the first element, as they're the resuming 
            # of every cohort
            AN = int(ANs[0][-1])
            AC = int(ACs[0][-1])
            try:
                AF = float(AC/AN)
            except ZeroDivisionError:
                AF = 0
            return AN, AC, AF
        vcf = cyvcf2.VCF(os.path.join(self.outfolder, self.main_filename), threads=4)
        vcf.add_info_to_header({'ID': 'AN', 'Description': 'Total allele in genotypes', 'Type':'Integer', 'Number': 'A'})
        vcf.add_info_to_header({'ID': 'AC', 'Description': 'Allele count in genotypes', 'Type':'Integer', 'Number': 'A'})
        vcf.add_info_to_header({'ID': 'AF', 'Description': 'Allele Frequency', 'Type':'Float', 'Number': 'A'})
        # open the writer
        fpath = os.path.join(self.outfolder, self.main_filename.replace(".gz", "_withAF.gz"))
        w = cyvcf2.Writer(fpath, vcf, mode="wz")
        for variant in vcf:
            variant.INFO['AN'], variant.INFO['AC'], variant.INFO['AF'] = calc_freq(variant)
            w.write_record(variant)
        vcf.close()
        w.close()
        # do also the index
        subprocess.run([
            'tabix', '-p', 'vcf', fpath
        ])
        # update the target filename accordingly
        self.main_filename = fpath.split('/')[-1]




def update_yaml(conf_main: str, resources: str, outfolder: str):
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
    for res_name in conf_main_yaml["resources"]:
        # first ensure that the file is not present. 
        if not os.path.isfile(conf_main_yaml["resources"][res_name]):
            if not res_name in resources.keys():
                print(f"Unable to find the URL for {res_name}")
                continue
            else:
                # download regularly 
                resource_entry = ResourceEntry(
                    conf_file=conf_main_yaml,
                    resources_entry=resources[res_name],
                    res_name=res_name,
                    outfolder=outfolder)
                # genome, transcriptome and gtf doesn't need conversions
                if resource_entry.res_type in ["fasta", "gtf"]:
                    resource_entry._download_stuff()
                elif resource_entry.res_type == "VCF":
                    # check first the notation, then proceed as needed
                    chr_notation = ChromosomeConverter.infer_notation_vcf(
                        vcf_url=resource_entry.resources_entry["url"]
                    )
                    if chr_notation == "ensembl":
                        # no conversion. just download
                        resource_entry._download_stuff()
                        # do the index
                        subprocess.run([
                            'tabix', '-p', 'vcf', os.path.join(outfolder, resource_entry.main_filename)
                        ])
                    else:
                        # convert on the fly with bcftools and the adequate file
                        cmd1 = f"bcftools view {resource_entry.resources_entry['url']} | bcftools annotate --rename-chrs {os.path.join(outfolder, chr_notation)}_to_ensembl.csv | bgzip -c > {os.path.join(outfolder, resource_entry.main_filename)}"
                        subprocess.run([cmd1], shell=True)
                        subprocess.run(["tabix", "-p", "vcf", os.path.join(outfolder, resource_entry.main_filename)])
                        # if dbSNPs, we need to adjust also for AF 
                        if resource_entry.res_name == "dbsnps":
                            resource_entry.generate_allele_frequency()
                elif resource_entry.res_type == "table":
                    # that's the stuff we need to do for the REDI portal file
                    chr_converter.convert_REDI(
                        bed_url=resource_entry.resources_entry["url"],
                        bed_output=os.path.join(outfolder, "REDI_portal.BED")
                    )
                    resource_entry.main_filename = "REDI_portal.BED.gz"
                elif resource_entry.res_type == "archive":
                    # that's for VEP: download and extract
                    resource_entry._download_stuff()
                    subprocess.run([
                        'tar', '-xzvf', os.path.join(outfolder, resource_entry.main_filename)
                    ])
                    resource_entry.main_filename = resource_entry.main_filename.replace('.tar.gz', '')
                # update entry accordingly
                conf_main_yaml[res_name] = os.path.join(outfolder, resource_entry.main_filename)
                
    # that's good, now we could write out the YAML 
    with open(conf_main, "w") as conf_main:
        yaml.dump(conf_main_yaml, conf_main)

if __name__ == "__main__":
    chr_converter = ChromosomeConverter()
    # generate table
    chr_converter.generate_txt("ucsc", "ensembl", "ucsc_to_ensembl.csv")
    update_yaml(conf_main=sys.argv[1], resources=sys.argv[2], outfolder=sys.argv[3])

