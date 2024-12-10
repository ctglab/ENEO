import cyvcf2
import pandas as pd
import subprocess
import os
import time
import glob
import argparse
import re
import numpy as np
import shutil
from multiprocessing import Pool


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', type=str, help='VCF file')
    parser.add_argument('-p', '--patient', type=str, help='Patient name')
    parser.add_argument('-o', '--outfolder', type=str, help='Output folder')
    parser.add_argument('-a', '--alleles', type=str, help='Single line allele file with alleles separated by commas')
    parser.add_argument('-n', '--ncpus', type=int, default=10, help='Number of CPUs to use')
    return parser.parse_args()

class variantCollector(object):
    def __init__(self, vcf_file: cyvcf2.VCF, patname) -> None:
        self.variants = self._store_variants(vcf_file)
        self.patname = patname 
    
    def _store_variants(self, vcf_file: cyvcf2.VCF, onlyAutosomicOrSexual=True):
        # the objective here is to store variants to associate also 
        # eventual phased variants.
        variants = {}
        for variant in vcf_file:
            if onlyAutosomicOrSexual:
                if str(variant.CHROM) in [f"chr{i}" for i in range(23)] + ['X', 'Y']:
                    variants[f'{variant.CHROM}:{variant.POS}'] = variant
            else:
                variants[f'{variant.CHROM}:{variant.POS}'] = variant
        return variants
                    
class variantExtended(object):
    def __init__(self, variant: cyvcf2.Variant):
        csq_keys = ["Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","TSL","FrameshiftSequence","WildtypeProtein"]
        self.variant = variant
        self.position = f"{self.variant.CHROM}:{self.variant.POS}-{self.variant.POS + 1}"
        self.germProb = 1 - self.variant.INFO.get('somProb')
        if self.variant.INFO.get('CSQ') is not None:
            self.csq = dict(zip(csq_keys, self.variant.INFO.get('CSQ').split('|')))
            try:
                self.aa_change, self.cdna_change = self._get_changes()
            except KeyError:
                print(self.csq)
                exit()
        else:
            self.aa_change, self.cdna_change = None, None
            self.csq = None
        self.phased_variants = self._get_phased()
        self.stringID = self._get_stringID()
        
    
    def _get_changes(self):
        if not self.csq['IMPACT'] in ['MODERATE', 'HIGH']:
            return None, None
        else:
            aa_change = f"p.{self.csq['Amino_acids'].split('/')[0]}{self.csq['Protein_position']}{self.csq['Amino_acids'].split('/')[1]}"
            cDNA_change = f"c.{self.variant.REF}{self.csq['cDNA_position']}{self.variant.ALT[0]}"
            return aa_change, cDNA_change
    def _get_phased(self):
        if self.variant.format('PS') is None:
            return []
        else:
            phased = []
            for variant in self.variant.format('PS')[0]:
                phased.append(f"{self.variant.CHROM}:{variant}")
            return phased
    
    def _get_substrings_at_index(self, s, index, length):
        substrings = []
        for start in range(max(0, index - length + 1), min(index + 1, len(s) - length + 1)):
            substrings.append(s[start:start + length])
        return substrings


    @staticmethod   
    def _get_sequences(self,  collector: variantCollector, min_length=8, max_length=12):
        wt_sequences = []
        mt_sequences = []
        aa_changes = []
        # control for the presence of any other phased variant
        if len(self.phased_variants) == 0:
            aa_changes.append(self.aa_change)
            mt_seq = self.csq["WildtypeProtein"][:int(self.csq['Protein_position']) - 1] + self.aa_change[-1] + self.csq["WildtypeProtein"][int(self.csq['Protein_position']):]
        else:
            # here we control for the presence of any other phased variant. We used the variant collector to 
            # get back infos about the other variant, then use the information about the modification on the 
            # aminoacidic sequence to generate the new mt_sequence which will be used for slicing.
            edits = {}
            # first add the current variant
            edits[self.csq['Protein_position']] = self.csq['Amino_acids'].split('/')[1]
            for phased_var in self.phased_variants:
                try:
                    phased_variant = collector.variants[phased_var]
                    phased_variant_extended = variantExtended(phased_variant)
                    # not all the variants have an impact on the aa sequence
                    if phased_variant_extended.aa_change is not None:
                        edits[phased_variant_extended.csq['Protein_position']] = phased_variant_extended.csq['Amino_acids'].split('/')[1]
                        # add aa_changes to the collector
                        aa_changes.append(phased_variant_extended.aa_change)
                    else:
                        continue
                except KeyError:
                    # this is the case where a phased variant was filtered
                    continue
            mt_seq = self.csq["WildtypeProtein"]
            for edit in edits:
                mt_seq = mt_seq[:int(edit) - 1] + edits[edit] + self.csq["WildtypeProtein"][int(edit):]
        #now generates frame within the sequence
        index = int(self.csq['Protein_position']) - 1
        for length in range(min_length, max_length + 1):
            wt_sequences.extend(self._get_substrings_at_index(self.csq['WildtypeProtein'], index, length))
            mt_sequences.extend(self._get_substrings_at_index(mt_seq, index, length))
        return wt_sequences, mt_sequences, aa_changes
            
    def _get_stringID(self):
        return f"{self.variant.CHROM}:{self.variant.POS}"

class frameGenerator(object):
    def __init__(self, hlas: list, collector: variantCollector, min_length=8, max_length=12):
        self.hlas = hlas
        self.collector = collector
        self.min_length = min_length
        self.max_length = max_length
        self.frame = self._generate_frame()
    
    def _generate_frame(self, germProbthreshold=0.5):
        frame = []
        for variant in self.collector.variants.values():
            variant_extended = variantExtended(variant)
            if variant_extended.aa_change is not None:
                wt_seqs, mt_seqs, aa_changes = variantExtended._get_sequences(variant_extended, self.collector)
                # concatenate as a single string the aa_change
                aa_change = ','.join(aa_changes)
                for hla in self.hlas:
                    for wt_seq, mt_seq in zip(wt_seqs, mt_seqs):
                        frame.append({'patient': self.collector.patname, 'position': variant_extended.position, 'germProb': variant_extended.germProb, 'WT_Epitope_Seq': wt_seq, 'MT_Epitope_Seq': mt_seq, 'HLA': hla, 'aa_change': aa_change })
        frame = pd.DataFrame(frame)
        # drop duplicates within the frame
        frame = frame.drop_duplicates(subset=['MT_Epitope_Seq', 'HLA'])
        # subset the frame
        frame = frame[frame['germProb'] <= germProbthreshold]
        return frame
    
def convertHLA_notation(hla: str):
        """
        netMHCpan doesn't like the HLA notation used by the other tools.
        """
        if '*' in hla:
            return hla.replace('*', '')
        else:
            return f"{hla[0:5]}*{hla[-5:]}"

class launcher(object):
    def __init__(self, inputDf: pd.DataFrame, seqColumn: str, alleleCol: str, outfolder: str, ncpus=10, **kwargs):
        if any([seqColumn not in inputDf.columns, alleleCol not in inputDf.columns]):
            raise ValueError("Column names not found in input dataframe")
        else:
            self.sequences = list(inputDf[seqColumn].values)
            self.alleles = list(inputDf[alleleCol].values)
            self.df = pd.DataFrame({'pep': self.sequences, 'mhc': self.alleles})
            if len(self.sequences) != len(self.alleles):
                raise ValueError("Number of sequences and alleles do not match. Check input file")
            self.outfolder = outfolder
            self.ncpus = ncpus
            self.params = kwargs
    
class netMHCpan_launcher(launcher):
    def __init__(self, inputDf: pd.DataFrame, seqColumn: str, alleleCol: str, outfolder:str,  **kwargs):
        super().__init__(inputDf, seqColumn, alleleCol, outfolder,**kwargs)
        self.isok = self.ensureConfig()
        self.exec_path = self.params['exec_path']
        self.cmd_params = self.parse_params()

    def run_subprocess(self, args):
        exec_path, seqfile, conv_allele, cmd_params, outfile = args
        with open(outfile, 'w') as ofile:
            subprocess.run([exec_path, '-p', seqfile, '-a', conv_allele, cmd_params], stdout=ofile)
    
    def ensureConfig(self):
        """
        Simple check to verify that the netmhcpan is correcly configured. 
        @TODO: just grab from the exec path
        """
        if not os.path.exists(self.params['data_path']):
            raise ValueError("netMHCpan data path not found. Please download them and place in the right folder.")
        else:
            return True
            
    def parse_params(self):
        """
        Using the external parameters passed through the class, return as a string
        to be used for the subprocess
        """
        string_params = ''
        if 'rth' in self.params.keys():
            string_params += f"-rth {self.params['rth']}"
        if 'rlt' in self.params.keys():
            string_params += f"-rlt {self.params['rlt']}"
        return string_params
    
    def launch(self):
        # Multiprocessing version
        if self.isok:
            print("Executing netmhcpan..")
            # create the output directory for the results
            start = time.time()
            outfolder = self.outfolder
            os.makedirs(outfolder, exist_ok=True)
            # create the input file
            # make also a temp folder 
            tempfolder = os.path.join(outfolder, 'temp')
            os.makedirs(tempfolder, exist_ok=True)
            # split unique_alleles into chunks
            n_chunks = self.ncpus
            chunks = np.array_split(self.df, n_chunks)
            columns = self.df.columns
            args_list = []
            for idx,chunk in enumerate(chunks):
                print(f"Processing chunk {idx} of {n_chunks}")
                
                chunk = pd.DataFrame(chunk)
                chunk.columns = columns
                for allele in chunk['mhc'].unique():
                    conv_allele = convertHLA_notation(allele)
                    seqs = chunk[chunk['mhc'] == allele]['pep']         
                    seqfile = os.path.join(tempfolder, f'{conv_allele}_batch_{idx}.csv')
                    seqs.to_csv(seqfile, index=False, header=False)
                    outfile = os.path.join(tempfolder, f'{conv_allele}_results_batch_{idx}.csv')
                    args_list.append((self.exec_path, seqfile, conv_allele, self.cmd_params, outfile))
            # launch
            with Pool() as pool:
                pool.map(self.run_subprocess, args_list)
            end = time.time()
            print(f"netMHCpan completed the task in {end-start:.3f} seconds.")

    def parse_output(self):
        """
        Parse the output of netMHCpan and return a dataframe
        """
        temp_dest_path = os.path.join(self.outfolder, 'temp')
        files_to_parse = glob.glob(temp_dest_path+'/*_results*.csv')
        dest_df = pd.DataFrame()
        # get alleles and batches all together
        alleles = list(set([x.split('/')[-1].split('_')[0] for x in files_to_parse]))
        for allele in alleles:
            conv_allele = convertHLA_notation(allele)
            # subset files
            subset_files = [x for x in files_to_parse if allele in x]
            outfile = os.path.join(self.outfolder, f'{allele}.csv')
            with open(outfile, 'w') as ofile:
                for f in subset_files:
                    # get the batch number
                    batch_i = f.split('/')[-1].split('_')[-1].split('.')[0]
                    original_input = os.path.join(self.outfolder, 'temp', f'{allele}_batch_{batch_i}.csv')
                    # get the number of sequences
                    num_seqs = 0
                    with open(original_input, 'r') as original_iput:
                        num_seqs = len(original_iput.readlines())
                    with open(f, 'r') as infile:
                        good_header = 0
                        bad_header = 0
                        for idx,line in enumerate(infile):
                            # remove whitespaces and replace with commas
                            edited_line = re.sub(r'\s+', ',', line.strip())
                            if edited_line.startswith('Pos'):
                                # good. this is the header.
                                # plus one because of the empty line
                                good_header = idx
                                bad_header = idx + num_seqs + 1
                                ofile.write(edited_line + '\n')
                            elif idx > good_header and idx <= bad_header:
                                if not line.startswith('-'):
                                    # some lines are offending for the parsing. Damn how on hell someone 
                                    # can write output file like this one?
                                    if len(edited_line.split(',')) == 15:
                                        # the last two field are just the same one, but with a space between.
                                        # just connect them into a single one
                                        last_field = ' '.join(edited_line.split(',')[-2:])
                                        edited_line = ','.join(edited_line.split(',')[0:-2] + [last_field])
                                    ofile.write(edited_line + '\n')
                            else:
                                continue
            try:
                df = pd.read_csv(outfile, skipinitialspace=True, index_col=False)
            except:
                print("Something went wrong with the allele {}".format(allele))
                exit()
            try:
                dest_df = pd.concat([dest_df, df], ignore_index=True)
                # then we could remove the file
                os.remove(outfile)
            except KeyError:
                print("Something went wrong with the allele {}".format(allele))
                exit() 
        #drop temp 
        shutil.rmtree(temp_dest_path)
        # os.remove(os.path.join(self.outfolder, 'temp'))
        dest_df = dest_df.loc[:,['HLA','Peptide','Score_EL','Rnk_EL']].rename(columns={
            'Peptide': 'MT_Epitope_Seq', 'Score_EL': 'score_EL_netmhcpan', 'Rnk_EL': 'rank_EL_netmhcpan'})   
        return dest_df


if __name__ == "__main__":
    args = parse_args()
    vcf = cyvcf2.VCF(args.vcf)
    hla = open(args.alleles, 'r').readlines()[0].strip().split(',')
    variants = variantCollector(vcf, patname=args.patient)
    print(len(variants.variants))
    frame = frameGenerator(hla, variants).frame
    # instantiate the netMHCpan launcher
    netMHCpan = netMHCpan_launcher(frame, 'MT_Epitope_Seq', 'HLA', args.outfolder, exec_path="/opt/netmhcpan/netMHCpan", data_path="/opt/netmhcpan/data", ncpus=args.ncpus)
    netMHCpan.launch()
    netMHCpan_df = netMHCpan.parse_output()
    output_frame = pd.merge(frame, netMHCpan_df, on=['MT_Epitope_Seq', 'HLA'], how='left')
    output_frame.to_csv(os.path.join(args.outfolder, f'{args.patient}.epitopes.csv'), index=False)