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

def phred_to_prob(phred):
    return 10 ** (-phred / 10)

def prob_to_phred(prob):
    return -10 * np.log10(prob)


class popFreq(object):
    def __init__(self, info_field):
        self.default_freq = 0.01/(71702*2)
        self.info_field = dict(info_field)
        self.info_field = {k.lower():v for k,v in self.info_field.items()}
        self.populations = set([x.split('_')[0] for x in self.info_field.keys() if all([x.endswith('_af'), "ex" not in x])])
        if len(self.populations) > 0:
            self.AFs = {pop: self.info_field[f'{pop}_af'] for pop in self.populations}
            self.ANs = {pop: self.info_field[f'{pop}_an'] for pop in self.populations}
            self.ACs = {pop: self.info_field[f'{pop}_ac'] for pop in self.populations}
            # fix the gnomad
            self._fix_info()
            self.avg_freq = self._get_avg_freq()
        else:
            self.avg_freq = self.default_freq
    
    def _fix_info(self):
        # this function has to handle a lot of scenarios. The relevant info 
        # is the AF here, so we control the size of other fields using the 
        # AF field.
        to_drop = []
        for pop in self.populations:
            if isinstance(self.AFs[pop], (int, float)):
                # we're in a single allelic variant, found in the population
                if isinstance(self.ANs[pop], (int, float)):
                    # nothing to change
                    pass
                elif isinstance(self.ANs[pop], tuple):
                    # more than a single entry here: just return the non empty one
                    try:
                        an = [x for x in self.ANs[pop] if x != None][0]
                        self.ANs[pop] = an
                    except IndexError:
                        # all the values for AN are empty
                        raise("All the values of AN are None. Something to check here")
                else:
                    #@TODO: we may want to update then using ACs
                    raise ValueError("AN could not be none when AF is computed.")
            elif isinstance(self.AFs[pop], tuple):
                # we're in a multiallelic variant. Note that the we expect to have a single value 
                # given the annotation logic.
                try:
                    af = [x for x in self.AFs[pop] if x != None][0]
                    self.AFs[pop] = af
                    # update also the ANs
                    if isinstance(self.ANs[pop], tuple):
                        self.ANs[pop] = [x for x in self.ANs[pop] if x!= None][0]
                except IndexError:
                    # all the values are none for AF. 
                    if isinstance(self.ANs[pop], (int, float)):
                        # this is a strange case. However, it's safe to reconduct this to 
                        # a 0 AF for this variant
                        self.AFs[pop] = 0
                    elif isinstance(self.ANs[pop], tuple):
                        # more than a AN, usually equal each other or one of them is None.
                        try:
                            self.ANs[pop] = [x for x in self.ANs[pop] if x != None][0]
                            self.AFs[pop] = 0
                        except IndexError:
                            # we have none for both. This is then a strange case: we're going to remove 
                            # the population from the dict
                            to_drop.append(pop)
                    else:
                        # the AN is None, but we'd multiple values for AFs
                        # this is again a strange case
                        #@TODO: should we try to recover this from ACs?
                        raise ValueError("How could you have an AF with None in AN????")
            else:
                # AF is none
                if isinstance(self.ANs[pop], (int, float)):
                    # just put Os for the AF here
                    self.AFs[pop] = 0
                elif isinstance(self.ANs[pop], tuple):
                    # multiple values for ANs, man
                    try:
                        self.ANs[pop] = [x for x in self.ANs[pop] if x != None][0]
                        self.AFs[pop] = 0
                    except IndexError:
                        # we have none for both. This is then a strange case: we're going to remove 
                        # the population from the dict
                        to_drop.append(pop)
                else:
                    to_drop.append(pop)
        if len(to_drop) > 0:
            for pop in to_drop:
                self.populations.remove(pop)
                del self.AFs[pop]
                del self.ANs[pop]

    def _get_avg_freq(self):
        # we adjust before the infos, so we'd just to return the average here
        ANs = [self.ANs[pop] for pop in self.populations]
        AFs = [self.AFs[pop] for pop in self.populations]
        try:
            return np.average(
                AFs, weights=ANs
            ) 
        except ZeroDivisionError:
            return self.default_freq 
        except TypeError:
            raise

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
        self.avg_freq = popFreq(self.variant.INFO).avg_freq
        self.germ_prob = self._germProb()
                
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
        post0_1 = np.round(p0_1 / (p_S + p0_1 + p1_1), 3)
        post1_1 = np.round(p1_1 / (p_S + p0_1 + p1_1), 3)
        p_S = np.round(p_S / (p_S + p0_1 + p1_1), 3)
        return (post0_1, post1_1, p_S)


def writeGermProb(input_vcf: str, output_vcf: str):
    """
    This function writes the germProb to the VCF file
    """
    vcf_file = cyvcf2.VCF(input_vcf, threads=8)
    vcf_file.add_info_to_header(
        {'ID': 'hetProb', 'Description': 'Probability germline for alt heterozygosity given population allele frequency', 'Type': 'Float', 'Number': "1"})
    vcf_file.add_info_to_header(
        {'ID': 'homaltProb', 'Description': 'Probability of germline for alt homozygosity given population allele frequency', 'Type': 'Float', 'Number': "1"})
    vcf_file.add_info_to_header(
        {'ID': 'somProb', 'Description': 'Probability of somatic event for all genotypes given population allele frequency', 'Type': 'Float', 'Number': "1"}
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