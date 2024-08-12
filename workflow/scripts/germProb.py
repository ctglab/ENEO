import cyvcf2
import numpy as np
import scipy.stats as stats
import sys

########### PRIORS ###########
# For concept behind this calculation, refer to
# the accompanying publication
##############################

alpha_exonMut, beta_exonMut = (0.87, 612.23)
alpha_MutRatio, beta_MutRatio = (0.95, 199.58)

##################################


def get_avg_freq(variant: cyvcf2.Variant):
    infodict = dict(variant.INFO)
    ANs = []
    AFs = []
    for AF, AN in zip(
        [x for x in infodict.keys() if 'AF' in x],
        [x for x in infodict.keys() if 'AN' in x]):
        ANs.append(infodict[AN])
        AFs.append(infodict[AF])
    # drop nan in both arrays
    ANs = [x for x in ANs if isinstance(x, float) or isinstance(x, int)]
    AFs = [x for x in AFs if isinstance(x, float) or isinstance(x, int)]
    if len(ANs) == len(AFs) and len(AFs) > 0:
        try:
            return np.average(AFs, weights=ANs)
        except ZeroDivisionError:
            return 0
    else:
        return np.nan

def phred_to_prob(phred: float):
    return 10**(-phred/10)


def prob_to_phred(prob: float):
    return -10*np.log10(prob)

class variantProb():
    """
    This class extends the cyvcf2.Variant class to include the probability of being a germ-line variant.
    """
    def __init__(self, variant: cyvcf2.Variant):
        self.variant = variant
        self.fixed_pi = stats.beta.rvs(alpha_exonMut, beta_exonMut) * stats.beta.rvs(alpha_MutRatio, beta_MutRatio)
        self.is_het = self._read_genotype()[0]
        self.is_homalt = self._read_genotype()[1]
        self.is_homref = self._read_genotype()[2]
        self.is_haploid = self._read_genotype()[3]
        self.is_phased = self._read_genotype()[4]
        self.hetLT = self._get_LTs("het")
        self.homaltLT = self._get_LTs("homalt")
        self.avg_freq = self._get_avg_freq()
        self.germ_prob = self._germProb()

    def _get_avg_freq(self):
        infodict = dict(self.variant.INFO)
        ANs = []
        AFs = []
        for AF, AN in zip(
            [x for x in infodict.keys() if 'AF' in x],
            [x for x in infodict.keys() if 'AN' in x]):
            ANs.append(infodict[AN])
            AFs.append(infodict[AF])
        # drop nan in both arrays
        ANs = [x for x in ANs if isinstance(x, float) or isinstance(x, int)]
        AFs = [x for x in AFs if isinstance(x, float) or isinstance(x, int)]
        if len(ANs) == len(AFs) and len(AFs) > 0:
            try:
                return np.average(AFs, weights=ANs)
            except ZeroDivisionError:
                return np.round(0.01/(71702*2), 4)
        else:
            # when we have no AF or AN, we use the priors as defined in the Mutect2 paper.
            # they stated to use the mean of the posterior defined by 
            # Beta(α, β + N), which is basically α/N
            # where N is equal to gnomAD sample size * 2 (chromosomes)
            return np.round(0.01/(71702*2), 4)
                
    def _read_genotype(self):
        """
        Just an helper to handle haploid scenario. Motivated by
        https://github.com/brentp/cyvcf2/issues/204
        """
        # unpack genotype
        is_het = False
        is_homalt = False
        is_homref = False
        is_haploid = False
        phased = False
        fields = self.variant.genotypes[0]
        if len(fields) == 2:
            # haploid call
            is_haploid = True
            is_homalt = True if fields[0] == 1 else False 
            is_homref = True if fields[0] == 0 else False
        else:
            # diploid call
            is_homalt = True if fields[0] == 1 and fields[1] == 1 else False
            is_het = True if fields[0] != fields[1] else False
            is_homref = True if fields[0] == 0 and fields[1] == 0 else False
        return (is_het, is_homalt, is_homref, is_haploid, phased)

    def _get_LTs(self, which="het") -> float:
        """
        Helper with the same concept as before

        Args:
            which (str): "het" or "homalt"
        """
        if self.is_haploid:
            # this completely break the assignation of objects for the cyvcf2 API.
            if self.is_homalt:
                if which == "het":
                    return prob_to_phred(np.finfo(float).eps)
                else:
                    return self.variant.gt_phred_ll_het[0]
            else:
                raise NotImplementedError("Haploid homref not implemented")
        else:
            if which == "het":
                return self.variant.gt_phred_ll_het[0]
            else:
                return self.variant.gt_phred_ll_homalt[0]


    def __str__(self):
        return str(self.variant)


    def _germProb(self) -> tuple:
        p_S = self.fixed_pi
        p0_1 = phred_to_prob(self.hetLT) * self.avg_freq * (1 - self.avg_freq) 
        p1_1 = self.avg_freq**2 * phred_to_prob(self.homaltLT) 
        # get normalized probabilities
        post0_1 = prob_to_phred(p0_1 / (p_S + p0_1 + p1_1))
        post1_1 = prob_to_phred(p1_1 / (p_S + p0_1 + p1_1))
        p_S = prob_to_phred(p_S / (p_S + p0_1 + p1_1))
        return (post0_1, post1_1, p_S)

def writeGermProb(input_vcf: str, output_vcf: str):
    """
    This function writes the germProb to the VCF file
    """
    vcf_file = cyvcf2.VCF(input_vcf, threads=8)
    vcf_file.add_info_to_header(
        {'ID': 'hetProb', 'Description': 'PHRED score of germline for alt heterozygosity given population allele frequency', 'Type': 'Float', 'Number': "1"})
    vcf_file.add_info_to_header(
        {'ID': 'homaltProb', 'Description': 'PHRED score of germline for alt homozygosity given population allele frequency', 'Type': 'Float', 'Number': "1"})
    vcf_file.add_info_to_header(
        {'ID': 'somProb', 'Description': 'PHRED score of somatic event for all genotypes given population allele frequency', 'Type': 'Float', 'Number': "1"}
    )
    outfile = cyvcf2.Writer(output_vcf, tmpl=vcf_file, mode='wz')
    for variant in vcf_file:
        variant_adj = variantProb(variant)
        variant.INFO['hetProb'] = variant_adj.germ_prob[0]
        variant.INFO['homaltProb'] = variant_adj.germ_prob[1]
        variant.INFO['somProb'] = variant_adj.germ_prob[2]
        outfile.write_record(variant)
    outfile.close()
    

if __name__ == "__main__":
    writeGermProb(sys.argv[1], sys.argv[2])