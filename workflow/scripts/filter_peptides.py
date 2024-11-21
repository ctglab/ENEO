import argparse
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description='Filter predictions')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output file')
    parser.add_argument('-c', '--calibration_data', type=str, required=True, help='Path to the calibration data')
    parser.add_argument('-a', '--hla_ligand_atlas', type=str, required=True, help='Path to the HLA ligand atlas')
    parser.add_argument('-min', '--min_length', type=int, default=8, help='Minimum length of the epitope')
    parser.add_argument('-max', '--max_length', type=int, default=12, help='Maximum length of the epitope')
    parser.add_argument('-g', '--germProb', type=float, default=.5, help='Maximum germProb')
    return parser.parse_args()

def filterPredictions(
    df: pd.DataFrame, 
    calibration_data: pd.DataFrame,
    HLA_ligand_atlas: pd.DataFrame,
    min_length: int=8,
    max_length: int=12,
    germProb: float=.5):
    """
    Several filtering steps are done here.
    """ 
    df.dropna(subset=['rank_EL_netmhcpan'], inplace=True)
    # filter by length and germProb
    df["len"] = df['MT_Epitope_Seq'].str.len()
    df = df.loc[(df['len'] >= min_length)&(df['len'] <= max_length)]
    df.drop(columns={'len'}, inplace=True)
    df = df.loc[df['germProb'] <= germProb]
    # filter by percentile
    df = df.merge(calibration_data.loc[:, ['HLA', 'optimal_percentile']], left_on='HLA', right_on='HLA', how='left')
    # if no percentile is given, take 2 as a general accepted value
    df['optimal_percentile'] = df['optimal_percentile'].fillna(2)
    df = df.loc[df['rank_EL_netmhcpan'] <= df['optimal_percentile']].drop(columns=['optimal_percentile'])
    # now sort, and take for the given mut-HLA only the 2 best returning entries
    df.sort_values(by=['rank_EL_netmhcpan'], ascending=True, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df['rankPos'] = df.groupby(['aa_change', 'HLA']).cumcount(1)
    df = df.loc[df['rankPos'].isin([0,1])]
    # drop peptides inside the HLA_ligand_atlas 
    healthy_pept = HLA_ligand_atlas.loc[HLA_ligand_atlas['hla_class'].isin(['HLA-I', 'HLA-I+II'])].loc[:, ['peptide_sequence']]
    healthy_pept['healthy'] = True
    df = df.merge(healthy_pept, left_on=['MT_Epitope_Seq'], right_on=['peptide_sequence'], how='left')
    df['healthy'].fillna(False, inplace=True)
    df = df.loc[df['healthy'] == False]
    df.drop(columns=['healthy', 'peptide_sequence', 'rankPos'], inplace=True)
    # sort again and return
    df = df.reset_index(drop=True).sort_values(by=['rank_EL_netmhcpan'])
    return df

def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    calibration_data = pd.read_csv(args.calibration_data)
    hla_ligand_atlas = pd.read_csv(args.hla_ligand_atlas, sep='\t', compression='gzip')
    df = filterPredictions(df, calibration_data, hla_ligand_atlas, args.min_length, args.max_length, args.germProb)
    df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()