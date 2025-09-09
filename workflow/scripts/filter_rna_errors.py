#/bin/python3

"""
This script is responsible for multiple steps of filtering using a given RAW vcf coming from the 
variant calling. In details, this is responsible for the flagging of:

- Variants present in Mutect2 PoN (done using bedtools intersect)
- Variants present in GIAB complex regions (done using bedtools intersect)
- Drop of variants not in regions affecting the downstream protein
- Drop of leakage errors 
- Drop of variants too close to splicing regions

"""

import os
import sys
import subprocess
import argparse
import cyvcf2
import bionumpy as bnp
import numpy as np
import shutil
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ])

def parse_args():
    parser = argparse.ArgumentParser(description='Filtering of the VCF file to account for sites that may induce systematic errors')
    parser.add_argument('--vcf', help='VCF file')
    parser.add_argument('--pon', help='Mutect2 Panel of Normals')
    parser.add_argument('--giab', help='GIAB complex regions')
    parser.add_argument('--gtf', help='reference GTF file')
    parser.add_argument('--ref', help='Reference genome')
    parser.add_argument('--out', help='Output VCF file')
    parser.add_argument('--pat', help='Patient ID')
    parser.add_argument('--dist', help='Distance from splicing site', default=3)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

class GTF_record(object):
    """
    A GTF record is the first building block of the parser.
    The attribute field is parsed resulting in a dict.

    Attributes
    ----------
    chromosome : str
        The chromosome that the GTF record belongs to.
    source : str
        The source of the GTF record.
    feature_type : str
        The type of feature that the GTF record represents.
    start : int
        The start position of the feature in the chromosome.
    end : int
        The end position of the feature in the chromosome.
    score : str
        The score of the GTF record.
    strand : str
        The strand that the GTF record belongs to.
    phase : str
        The phase of the GTF record.
    length : int
        The length of the feature.
    attributes : dict
        A dictionary of attributes from the GTF record.

    Methods
    -------
    __init__(chromosome, source, feature_type, start, end, score, strand, phase, attributes)
        Initialize a GTF record object.
    parse_attributes(attributes)
        Parse the attributes of a GTF record.
    is_coding(feat_dict)
        Check if a GTF record is coding.
    """

    def __init__(
        self,
        chromosome,
        source,
        feature_type,
        start,
        end,
        score,
        strand,
        phase,
        attributes,
    ):
        self.chromosome = str(chromosome)
        self.source = source
        self.feature_type = feature_type
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.length = abs(self.end - self.start)
        self.attributes = GTF_record.parse_attributes(attributes)

    @staticmethod
    def parse_attributes(attributes):
        if isinstance(attributes, dict):
            return attributes
        else:
            attributes = attributes.replace('"', '').replace('; ', ';').replace(" ", "=")
            try:
                return {k:v for k,v in [x.split("=")[:2] for x in attributes.split(";")[:-1]]}
            except ValueError:
                print(attributes)
                raise


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
        try:
            if feat_dict["gene_biotype"] == "protein_coding":
                return True
            else:
                return False
        except KeyError:
            # that's something for ncRNA 
            return False

class SplicingCollector(object):
    def __init__(self, gtf_file: str):
        self.gtf_file = gtf_file
        self.ss = self._get_splicing_sites()
    
    def _get_splicing_sites(self):
        ss = {}
        with open(self.gtf_file, "r") as gtf:
            for line in gtf:
                if not line.startswith("#"):
                    record = GTF_record(*line.rstrip().split("\t"))
                    if not record.chromosome in ss:
                        ss[record.chromosome] = set()
                    if record.feature_type == "exon" and GTF_record.is_coding(record.attributes):
                        #@TODO: we may want to retain just splicing junctions of transcript with TSL=1
                        # as this may interfere with the junctions of other transcripts that are less supported.
                        if record.attributes["exon_number"] == 1:
                            ss[record.chromosome].add(record.start)
                        else:
                            ss[record.chromosome].add(record.start)
                            ss[record.chromosome].add(record.end)
        # turn all the ss as a numpy array
        for k in ss:
            ss[k] = np.array(list(ss[k]))
        return ss

    def distance_from_splice(self, chrom: str, pos: int):
        if chrom not in self.ss:
            raise ValueError(f"Chromosome {chrom} not found in the splicing sites. Did you use the same chromosome notation?")
        else:
            return np.min(np.abs(self.ss[chrom] - pos))

class VariantCollector(object):
    def __init__(self, vcf_file: str, tag: str):
        self.vcf_file = vcf_file
        self.tag = tag
        self.variants = self._get_variants()
    
    def _get_variants(self):
        variants = []
        for variant in cyvcf2.VCF(self.vcf_file):
            variants.append(f"{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}")
        return set(variants)

class ProteinCodingCollector(object):
    def __init__(self, gtf: str):
        self.gtf = gtf
        self.regions = self._collect_regions()
    
    def _collect_regions(self):
        from collections import defaultdict
        intervals_dict = defaultdict(list)
        
        with open(self.gtf, "r") as gtf:
            for line in gtf:
                line = line.strip()
                if not line.startswith("#"):
                    record = GTF_record(*line.split("\t"))
                    # Only gather exons for coding transcripts with TSL=1
                    if record.feature_type == "exon" and GTF_record.is_coding(record.attributes):
                        if record.attributes.get("transcript_support_level") == "1":
                            intervals_dict[record.chromosome].append([record.start, record.end])

        # Merge and convert to sorted NumPy arrays for each chromosome
        merged = {}
        for chrom, intervals in intervals_dict.items():
            # Sort by start
            intervals.sort(key=lambda x: x[0])
            # Optionally merge overlapping intervals
            merged_intervals = []
            current_start, current_end = intervals[0]
            for i in range(1, len(intervals)):
                start, end = intervals[i]
                if start <= current_end:
                    # Overlapping, so extend
                    current_end = max(current_end, end)
                else:
                    merged_intervals.append([current_start, current_end])
                    current_start, current_end = start, end
            merged_intervals.append([current_start, current_end])

            # Convert to NumPy array
            merged[chrom] = np.array(merged_intervals, dtype=np.int64)
        return merged

    def is_in_region(self, chrom: str, pos: int):
        """
        Binary-search-based lookup.
        """
        if chrom not in self.regions:
            return False

        intervals_arr = self.regions[chrom]
        # intervals_arr is shape (N,2) with [start, end] sorted by start

        # Find the interval index whose start is just to the left of 'pos'
        idx = np.searchsorted(intervals_arr[:, 0], pos, side='right') - 1
        if idx < 0:
            return False

        # Check if pos <= intervals_arr[idx, 1]
        return pos <= intervals_arr[idx, 1]


def is_leakage_error(ref_string: str, alt: str):
    if alt == '-':
        return False
    if not len(ref_string) == 5:
        print(f"The reference {ref_string} is not 5 nucleotides long")
        return False
    try:
        if all([j==alt for j in ref_string[:2]]) and all([j==alt for j in ref_string[-2:]]):
            return True
        else:
            return False
    except ValueError:
        if len(alt) > 1:
            return False
        else:
            print(f"Error in the sequence {ref_string}")
            raise

def get_leakage_errors(vcf_file: str, genome: bnp.Genome):
    refs = []
    alts = []
    chroms = []
    starts = []
    ends = []
    positions = []
    leakage_errors = []
    errors_found = 0
    print("=="*20)
    print("Filtering leakage errors")
    for variant in cyvcf2.VCF(vcf_file):
        if variant.is_snp:
            refs.append(variant.REF)
            alts.append(variant.ALT[0])
            # get the interval from the variant
            chroms.append(variant.CHROM)
            starts.append(variant.POS - 3)
            ends.append(variant.POS + 2)
            positions.append(variant.POS)
    intervals = bnp.datatypes.Interval(chroms, starts, ends)
    # now we have the intervals, we can extract the sequences
    sequences = genome.read_sequence()[intervals]
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    # check now for leakage errors
    for i, seq in enumerate(sequences):
        if is_leakage_error(seq, alts[i]):
            leakage_errors.append(f"{chroms[i]}:{positions[i]}_{refs[i]}>{alts[i]}")
            errors_found += 1
    print(f"Total leakage errors: {errors_found}")
    return set(leakage_errors)

def get_variant_close_to_splice(vcf_file: str, sc: SplicingCollector, dist=4):
    print("=="*20)
    print(f"Filtering variants too close to splicing sites with a distance of {dist}")
    variants = []
    splicing_errors = 0
    for variant in cyvcf2.VCF(vcf_file):
        chrom = variant.CHROM
        pos = variant.POS
        # get the distance 
        distance = sc.distance_from_splice(chrom, pos)
        if distance < dist:
            variants.append(f"{chrom}:{pos}_{variant.REF}>{variant.ALT[0]}")
            splicing_errors += 1
    print(f"Total splicing errors: {splicing_errors}")
    return set(variants)


def is_filtered(variant: cyvcf2.Variant, collectors: list):
    for collector in collectors:
        tag = collector.tag
        if len(set.intersection(
            set([f"{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}"]),
            collector.variants
        )) > 0:
            return tag
    return None


def main():
    args = parse_args()
    # read the genome
    genome = bnp.Genome.from_file(args.ref)
    sc = SplicingCollector(args.gtf)
    # get the protein coding regions
    pc = ProteinCodingCollector(args.gtf)
    patient = args.pat
    # drop a temp folder
    tmp_folder = os.path.join(
        os.path.dirname(os.path.abspath(args.out)), args.pat)
    os.makedirs(tmp_folder, exist_ok=True)
    # intersect with the PoN
    subprocess.run(f"bedtools intersect -a {args.vcf} -b {args.pon} -header 2>/dev/null | bgzip -c > {tmp_folder}/{patient}_pon.vcf.gz", shell=True)
    subprocess.run(f"tabix -p vcf {tmp_folder}/{patient}_pon.vcf.gz", shell=True)
    PoN_variants = VariantCollector(f"{tmp_folder}/{patient}_pon.vcf.gz", "PoN")
    # intersect with the GIAB complex regions
    subprocess.run(f"bedtools intersect -a {tmp_folder}/{patient}_pon.vcf.gz -b {args.giab} -header 2>/dev/null | bgzip -c > {tmp_folder}/{patient}_giab.vcf.gz", shell=True)
    subprocess.run(f"tabix -p vcf {tmp_folder}/{patient}_giab.vcf.gz", shell=True)
    GIAB_variants = VariantCollector(f"{tmp_folder}/{patient}_giab.vcf.gz", "GIAB")
    # drop the variants not in the protein coding regions
    subprocess.run(f"bcftools view {args.vcf} -e 'INFO/CSQ != \"\"' -Oz -o {tmp_folder}/{patient}_nocsq.vcf.gz", shell=True)
    subprocess.run(f"tabix -p vcf {tmp_folder}/{patient}_nocsq.vcf.gz", shell=True)
    CSQ_variants = VariantCollector(f"{tmp_folder}/{patient}_nocsq.vcf.gz", "non-coding")
    # drop the leakage errors
    leakage_errors = get_leakage_errors(args.vcf, genome)
    print(sc.ss)
    splicing_errors = get_variant_close_to_splice(args.vcf, sc, dist=args.dist)
    # final writing
    vcf_in = cyvcf2.VCF(args.vcf)
    # we want to write this as a filter
    vcf_in.add_filter_to_header(
        {'ID': 'PoN', 'Description': 'Variant present in PoN', 'Type': 'Flag', 'Number': "1"})
    vcf_in.add_filter_to_header(
        {'ID': 'GIAB', 'Description': 'Variant present in GIAB complex regions', 'Type': 'Flag', 'Number': "1"})
    vcf_in.add_filter_to_header(
        {'ID': 'non-coding', 'Description': 'Variant not in protein coding region', 'Type': 'Flag', 'Number': "1"})
    vcf_in.add_filter_to_header(
        {'ID': 'leakage', 'Description': 'Variant is a leakage error', 'Type': 'Flag', 'Number': "1"})
    vcf_in.add_filter_to_header(
        {'ID': 'low-TLS', 'Description': 'Variant is not in a region with TLS=1', 'Type': 'Flag', 'Number': "1"})
    vcf_in.add_filter_to_header(
        {'ID': 'splicing', 'Description': f'Variant is too close to splicing site, using distance {args.dist}', 'Type': 'Flag', 'Number': "1"})
    out_vcf = cyvcf2.Writer(args.out, vcf_in)
    for variant in vcf_in:
        filter = is_filtered(variant, [PoN_variants, GIAB_variants, CSQ_variants])
        if filter is not None:
            variant.FILTER = filter
        else:
            if len(set.intersection(
                set([f"{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}"]),
                leakage_errors
            )) > 0:
                variant.FILTER = "leakage"
            elif len(set.intersection(
                set([f"{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}"]),
                splicing_errors
            )) > 0:
                variant.FILTER = "splicing"
            elif not pc.is_in_region(variant.CHROM, variant.POS):
                variant.FILTER = "low-TLS"
            else:
                variant.FILTER = "PASS"
        out_vcf.write_record(variant)
    out_vcf.close()
    # remove the temp folder
    shutil.rmtree(tmp_folder)
        

if __name__ == '__main__':
    main()
    sys.exit(0)