#!/usr/bin/python3 

import sys
import pandas as pd

def main(hla_genotype: str):
    """
    Simple function to read t1k output and return HLA types
    """
    df = pd.read_csv(hla_genotype, sep='\t', index_col=False, names=["gene_name","num_diff_alleles","allele_1","abundance_1","quality_1","allele_2","abundance_2","quality_2","secondary_alleles"])
    # filter to retain only HLA-A, HLA-B, HLA-C
    df = df[df["gene_name"].isin(["HLA-A","HLA-B","HLA-C"])]
    HLAs = []
    for row in df.itertuples():
        if row.quality_1 > 5:
            # if more than 1 allele is present, keep only the first
            if len(row.allele_1.split(",")) > 1:
                hla = ':'.join(row.allele_1.split(",")[0].split(":")[:-1])
            else:
                hla = ':'.join(row.allele_1.split(":")[:-1])
            if ":" in hla:
                HLAs.append(hla)

        if row.quality_2 > 5:
            if len(row.allele_2.split(",")) > 1:
                hla = ':'.join(row.allele_2.split(",")[0].split(":")[:-1])
            else:
                hla = ':'.join(row.allele_2.split(":")[:-1])
            if ":" in hla:
                HLAs.append(hla)
    print(','.join([x for x in set(HLAs)]))



if __name__ == '__main__':
    main(sys.argv[1])