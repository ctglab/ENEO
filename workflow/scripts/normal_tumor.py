#!/usr/bin/python3

import glob, sys
import pandas as pd

# Read the metadata file and return a Dataframe with
# sample	phenotype	patient


def samples_couple(metadata, patient):
    df = pd.read_csv(metadata)
    df["subject_id"] = df["subject_id"].map(lambda x: x.rsplit("NC"))
    df["phenotype"] = df["phenotype"].fillna("tumor")


df = pd.read_csv(sys.argv[1])
# transform the dataframe:
# leave just patient
df["subject_id"] = df["subject_id"].map(lambda x: x.rsplit("NC"))
df["phenotype"] = df["phenotype"].fillna("tumor")
