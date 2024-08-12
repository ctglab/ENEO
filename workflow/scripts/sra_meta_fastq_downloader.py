#!/usr/bin/python3

import glob
import sys
import os
from ftplib import FTP
import subprocess
import json


def get_json(project, accepted_projects):
    """
    Run ffq using provided project
    """
    if project[0:3] not in accepted_projects:
        print("Invalid project format provided.")
        print("Exiting..")
        sys.exit()
    else:
        try:
            subprocess.run(["ffq", "-t", project[0:3], project, "-o", ".", "--split"])
        except FileNotFoundError:
            print("ffq is not installed or may not be in $PATH" + "\n" + "Exiting..")
            sys.exit()


def get_metadata(project):
    """
    Create metadata from json, return list of ftp for fastq
    """
    with open(project + ".json", "r") as json_mdata:
        project_data = json.load(json_mdata)
        # this is a cool trick introduced in Python 3.5
        # Instead of doing `list(dict.keys())` it's possible to unpack into list any
        # iterable elements, like dict keys, using just `*`.
        accession_list = [*project_data["runs"]]
        ftp_url_list = []
        with open(project + "_metadata.csv", "w") as oput_file:
            oput_file.write("Sample,gender,phenotype,subject_id" + "\n")
            for x in accession_list:
                # json is a dict-nested structure
                sample_attr = project_data["runs"][x]["sample"]["attributes"]
                oput_file.write(
                    ",".join(
                        [
                            x,
                            sample_attr["gender"],
                            sample_attr["phenotype"],
                            str(sample_attr["subject_id"] + "\n"),
                        ]
                    )
                )
                seq_files = project_data["runs"][x]["files"]
                # if paired end, it's a list with length > 1
                if isinstance(seq_files, list) and len(seq_files) > 1:
                    # if paired, files is a list
                    for k in seq_files:
                        ftp_url_list.append(k["url"])
                else:
                    print("Something wrong with record " + x)
                    sys.exit()
    return ftp_url_list


def download_files(ftp_url_list):
    for x in ftp_url_list:
        print("Downloading file {}".format(x.split("/")[-1]))
        subprocess.run(["wget", "-c", x, "-P", "FASTQ"])
    print("Download finished.")
    print("Quitting..")
    sys.exit()


if __name__ == "__main__":
    accepted_project = ["SRR", "ERR", "DRR", "SRP", "ERP", "DRP", "GSE", "DOI"]
    try:
        project = sys.argv[1]
    except IndexError:
        print("No project provided. Exiting..")
        sys.exit()
    get_json(sys.argv[1], accepted_project)
    ftp_url_list = get_metadata(sys.argv[1])
    if len(ftp_url_list) > 1:
        download_choice = input(
            "I've found {} files. Would you like to download them? (Y/N)\t".format(
                len(ftp_url_list)
            )
        )
        if download_choice not in ["Y", "N"]:
            print("Only `Y` (yes) or `N` (no) are valid answers.")
            print("No is the default. Exiting..")
        elif download_choice == "Y":
            download_files(ftp_url_list)
        else:
            print("")
            print("Haven't you heard about my miserable past,")
            print("Critical past,")
            print("Pitiful past..")
            print("")
            sys.exit()
