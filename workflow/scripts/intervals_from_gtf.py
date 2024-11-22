import sys 
import requests

class GTF_record(object):
    def __init__(self,chromosome,source,feature_type,start,end,score,strand,phase,attributes):
        self.chromosome=str(chromosome)
        self.source=source
        self.feature_type=feature_type
        self.start=int(start) 
        self.end=int(end) 
        self.score=score
        self.strand=strand
        self.phase=phase
        self.length=abs(self.end - self.start)
        self.attributes = GTF_record.parse_attributes(attributes)

    @staticmethod    
    def parse_attributes(attributes):
        if isinstance(attributes, dict):
            return attributes
        else:
            feat_dict = {}
            keyVal = attributes.split(';')[:-1]
            keyVal = [x.lstrip().replace('"', '').replace(' ', '=') for x in keyVal]
            feat_dict = {k:v for k,v in [item.split('=') for item in keyVal]} 
            return feat_dict
        
    @staticmethod
    def is_coding(feat_dict: dict) -> bool:
        """
        Simple check if the given record is protein coding.

        Parameters
        ----------
        feat_dict : dict
            A dictionary of attributes from the GTF record.
        
        Returns
        -------
        bool
            True if the record is protein coding, False otherwise.
        """
        if feat_dict["gene_biotype"] == 'protein_coding':
            return True
        else:
            return False

    def asStr(self):
        """
        Return the record as a string for writing out.

        Returns
        -------
        str
            The record as a string, separated by tabs.
        """
        attributesAsStr = "; ".join([f'{k} "{v}"' for k,v in self.attributes.items()])
        return '\t'.join([str(x) for x in [self.chromosome, self.source, self.feature_type, self.start, self.end, self.score, self.strand, self.phase, attributesAsStr]]) + '\n'

def get_genes_from_kegg(pathway: str) -> list:
    """
    Get the genes from a KEGG pathway.

    Parameters
    ----------
    pathway : str
        The KEGG pathway to get the genes from.

    Returns
    -------
    list
        A list of genes from the KEGG pathway.
    """
    url = f"http://rest.kegg.jp/link/hsa/{pathway}"
    response = requests.get(url)
    if response.status_code == 200:
        genes = [x.split('\t')[1] for x in response.text.split('\n') if x]
        #now convert to gene symbols
        gene_symbols = []
        for gene in genes:
            url = f"https://rest.kegg.jp/get/{gene}"
            response = requests.get(url)
            # read the response
            if response.status_code == 200:
                gene_info = response.text.split('\n')
                for line in gene_info:
                    if line.startswith('SYMBOL'):
                        try:
                            gene_symbols.append(line.split('      ')[1].split(',')[0])
                            break
                        except IndexError:
                            print(line)
                            exit()
            else:
                print(f"Failed to get gene info for {gene}")
                exit()
        return gene_symbols
    else:
        raise ValueError(f"Failed to get genes from KEGG pathway {pathway}")

def main(gtf_file: str, outfile: str):
    pathway = "hsa04612"
    genes_to_exclude = get_genes_from_kegg(pathway)
    allowed_chrs = ['chrX', 'chrY'] + [f'chr{(i)}' for i in range(1, 23)]
    # keep also a stored list of sets to keep track of added intervals.
    intervals = set()
    with open(gtf_file, 'r') as gtf:
        with open(outfile, 'w') as outfile:
            for line in gtf:
                if line.startswith('#'):
                    continue
                else:
                    entry = GTF_record(*line.rstrip().split('\t'))
                    if entry.chromosome not in allowed_chrs:
                        continue
                    else:
                        if entry.is_coding and entry.feature_type == 'exon':
                            # we need also to drop too short exons
                            if entry.length < 10:
                                continue
                            else:
                                try:
                                    if entry.attributes['gene_name'] not in genes_to_exclude:
                                        intervals.add((entry.chromosome, entry.start, entry.end, entry.strand))
                                except KeyError:
                                    continue
            # cool, write out
            for interval in intervals:
                outfile.write('\t'.join([str(x) for x in interval]) + '\n')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
