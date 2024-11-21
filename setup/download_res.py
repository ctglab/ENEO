#!/usr/bin/python3

# This is an update of the legacy downloader, where we're moving into using the
# Gencode notation instead of the Ensembl one, to reduce the amount of conversions
# needed throughout the process.

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
        original_conv = pd.read_csv(
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz",
            compression="gzip",
            sep="\t",
            header=None,
            names=["source", "target", "source_type"],
        )
        return original_conv

    @staticmethod
    def infer_notation_vcf(self, vcf_url: str):
        cmd0 = f'bcftools view {vcf_url} | grep -v "#" | head -n 1 | cut -f1'
        chromosome_notation = subprocess.run(
            [cmd0], shell=True, stdout=subprocess.PIPE, encoding="utf-8"
        ).stdout.strip()
        try:
            inferred_not = self.conv_table.loc[
                self.conv_table["source"] == str(chromosome_notation)
            ]["source_type"].unique()[0]
            return inferred_not
        except IndexError:
            if chromosome_notation == "chr1":
                return "ucsc"
            else:
                raise NotImplementedError(
                    f"The chromosome {chromosome_notation} in unknown, and only homo sapiens is supported."
                )

    def generate_conv_file(self, chr_from: str, outfile=None):
        """
        Generate a simple TSV for chromosomes to move from a source notation to the Genbank one.

        Arguments:
        ------------
        chr_from: str
            Chromosome notation to move from. Accepted values are ["assembly", "ensembl", "genbank", "refseq"]

        """
        if chr_from.lower() not in ["ensembl", "assembly", "genbank", "refseq"]:
            raise ValueError(
                f"{chr_from} must be a choice between {', '.join(['ensembl', 'assembly', 'genbank', 'refseq'])}"
            )
        if outfile is None:
            outfile = os.path.join(os.path.dirname(cwd), f"{chr_from}_to_gencode.tsv")
        subset = self.conv_table.loc[
            self.conv_table["source_type"].str.contains(chr_from.lower())
        ]
        subset.loc[:, ["source", "target"]].to_csv(
            outfile, sep="\t", index=False, header=False
        )

    def convert_REDI(self, bed_url: str, bed_output: str, drop_intermediate=True):
        """
        Contrarly to the legacy method, this will just turn the file as a regular BED
        to be used for annotation
        """
        redi_df = pd.read_csv(
            bed_url, compression="gzip", sep="\t", usecols=[i for i in range(0, 9)]
        )
        # make it as a true BED, so converting with both start-end
        redi_df["Start"] = redi_df["Position"] - 1
        redi_df["End"] = redi_df["Start"] + 1
        redi_df = redi_df.loc[
            :,
            ["Region", "Start", "End"]
            + [
                x
                for x in redi_df.columns
                if x not in ["Region", "Position", "Start", "End"]
            ],
        ]
        redi_df.loc[
            :, ["Region"] + [x for x in redi_df.columns if x != "Region"]
        ].to_csv(bed_output, index=False, header=False, sep="\t")
        # sort & compress + index
        cmd0 = f"sort -k1,1 -k2,2n -k3,3n {bed_output} | bgzip -c > {bed_output+'.gz'}"
        subprocess.run([cmd0], shell=True)
        subprocess.run(["tabix", "-p", "bed", bed_output + ".gz"])
        if drop_intermediate:
            os.remove(bed_output)


class ResourceEntry(object):
    """
    This object handles the reading of the config file, check if file is present or
    has to be downloaded. If a conversion has to be done, it will be performed.
    Then it updates the relative configuration file accordingly.
    """

    def __init__(
        self, conf_file: dict, resources_entry: dict, res_name: str, outfolder: str
    ):
        self.conf_file = conf_file
        self.resources_entry = resources_entry
        self.res_name = res_name
        self.main_filename = self._get_main_filename()
        self.outfolder = outfolder
        self.downloaded = self._is_downloaded()
        self.res_type = self._get_restype()
        assert res_name in conf_file["resources"].keys()

    def _get_main_filename(self):
        return self.resources_entry["url"].split("/")[-1]

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
            subprocess.run(
                ["wget", "-c", self.resources_entry["url"], "-P", self.outfolder]
            )
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
            ANs = variant.format(field="AN")
            ACs = variant.format(field="AC")
            # the AN and AC fields are lists of strings
            # I'll take the first element, as they're the resuming
            # of every cohort
            AN = int(ANs[0][-1])
            AC = int(ACs[0][-1])
            try:
                AF = float(AC / AN)
            except ZeroDivisionError:
                AF = 0
            return AN, AC, AF

        vcf = cyvcf2.VCF(os.path.join(self.outfolder, self.main_filename), threads=4)
        vcf.add_info_to_header(
            {
                "ID": "AN",
                "Description": "Total allele in genotypes",
                "Type": "Integer",
                "Number": "A",
            }
        )
        vcf.add_info_to_header(
            {
                "ID": "AC",
                "Description": "Allele count in genotypes",
                "Type": "Integer",
                "Number": "A",
            }
        )
        vcf.add_info_to_header(
            {
                "ID": "AF",
                "Description": "Allele Frequency",
                "Type": "Float",
                "Number": "A",
            }
        )
        # open the writer
        fpath = os.path.join(
            self.outfolder, self.main_filename.replace(".gz", "_withAF.gz")
        )
        w = cyvcf2.Writer(fpath, vcf, mode="wz")
        for variant in vcf:
            variant.INFO["AN"], variant.INFO["AC"], variant.INFO["AF"] = calc_freq(
                variant
            )
            w.write_record(variant)
        vcf.close()
        w.close()
        # do also the index
        subprocess.run(["tabix", "-p", "vcf", fpath])
        # update the target filename accordingly
        self.main_filename = fpath.split("/")[-1]


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
                    outfolder=outfolder,
                )
                # genome, transcriptome and gtf doesn't need conversions
                if resource_entry.res_type in ["fasta", "gtf"]:
                    resource_entry._download_stuff()
                    # if it's a genome, do also the dictionary
                    if "genome" in resource_entry.main_filename:
                        outfile = os.path.join(outfolder, resource_entry.main_filename.replace('.fa.gz', '.dict'))
                        subprocess.run(
                            [
                                'gatk',
                                'CreateSequenceDictionary',
                                '-R',
                                os.path.join(outfolder, resource_entry.main_filename),
                                '-O',
                                outfile
                            ])
                elif resource_entry.res_name == "dbsnps":
                    # for dbsnps we need to do conversion and then annotation for allele frequency
                    refseq_conv_table = os.path.join(
                        os.path.dirname(cwd), "refseq_dbsnp.tsv"
                    )
                    chr_converter.generate_conv_file("refseq", refseq_conv_table)
                    cmd1 = f"bcftools annotate --rename-chrs {refseq_conv_table} {resource_entry.resources_entry['url']} | bgzip -c > {os.path.join(outfolder, resource_entry.main_filename)}"
                    subprocess.run([cmd1], shell=True)
                    subprocess.run(
                        [
                            "tabix",
                            "-p",
                            "vcf",
                            os.path.join(outfolder, resource_entry.main_filename),
                        ]
                    )
                    # now we need to add the allele frequency
                    resource_entry.generate_allele_frequency()
                elif resource_entry.res_type == "table":
                    # that's the stuff we need to do for the REDI portal file
                    chr_converter.convert_REDI(
                        bed_url=resource_entry.resources_entry["url"],
                        bed_output=os.path.join(outfolder, "REDI_portal.BED"),
                    )
                    resource_entry.main_filename = "REDI_portal.BED.gz"
                elif resource_entry.res_type == "archive":
                    # that's for VEP: download and extract
                    resource_entry._download_stuff()
                    subprocess.run(
                        [
                            "tar",
                            "-xzvf",
                            os.path.join(outfolder, resource_entry.main_filename),
                            "-C",
                            os.path.abspath(outfolder)
                        ]
                    )
                    resource_entry.main_filename = os.path.join(
                        outfolder, 'homo_sapiens'
                    )
                # update entry accordingly
                conf_main_yaml["resources"][res_name] = os.path.join(
                    os.path.abspath(outfolder), resource_entry.main_filename
                )

    # that's good, now we could write out the YAML
    with open(conf_main, "w") as conf_main:
        yaml.dump(conf_main_yaml, conf_main)


if __name__ == "__main__":
    chr_converter = ChromosomeConverter()
    # generate table
    update_yaml(conf_main=sys.argv[1], resources=sys.argv[2], outfolder=sys.argv[3])
