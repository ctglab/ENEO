# Calling intervals

By default, ENEO performs variant calling on exons of protein coding genes, with the exception of two types of calling regions:

- Known hard-to-call regions of the human genome. These regions were [extensively profiled](https://www.biorxiv.org/content/10.1101/2023.10.27.563846v1.full) by the GIAB consortium, and are [publicly available](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/). 

- Genes involved in the antigen presentation: given the initial goal of the pipeline to be applied in a personalized cancer vaccine setup, we discared these genes as a source of unwanted variations. We obtained the set of genes using the KEGG pathway annotated as 

??? example "Expand to see the full set of excluded genes in ENEO"
    | symbol   | description                                                                 |
    |----------|-----------------------------------------------------------------------------|
    | IFNG     | interferon gamma                                                            |
    | TNF      | tumor necrosis factor                                                       |
    | PSME1    | proteasome activator subunit 1                                              |
    | PSME2    | proteasome activator subunit 2                                              |
    | PSME3    | proteasome activator subunit 3                                              |
    | HSPA8    | heat shock protein family A (Hsp70) member 8                                |
    | HSPA1A   | heat shock protein family A (Hsp70) member 1A                               |
    | HSPA1L   | heat shock protein family A (Hsp70) member 1 like                           |
    | HSPA1B   | heat shock protein family A (Hsp70) member 1B                               |
    | HSPA6    | heat shock protein family A (Hsp70) member 6                                |
    | HSPA2    | heat shock protein family A (Hsp70) member 2                                |
    | HSPA4    | heat shock protein family A (Hsp70) member 4                                |
    | HSP90AA1 | heat shock protein 90 alpha family class A member 1                         |
    | HSP90AB1 | heat shock protein 90 alpha family class B member 1                         |
    | HLA-A    | major histocompatibility complex, class I, A                                |
    | HLA-B    | major histocompatibility complex, class I, B                                |
    | HLA-C    | major histocompatibility complex, class I, C                                |
    | HLA-F    | major histocompatibility complex, class I, F                                |
    | HLA-G    | major histocompatibility complex, class I, G                                |
    | HLA-E    | major histocompatibility complex, class I, E                                |
    | HSPA5    | heat shock protein family A (Hsp70) member 5                                |
    | CANX     | calnexin                                                                    |
    | B2M      | beta-2-microglobulin                                                        |
    | PDIA3    | protein disulfide isomerase family A member 3                               |
    | CALR     | calreticulin                                                                |
    | TAPBP    | TAP binding protein                                                         |
    | TAP1     | transporter 1, ATP binding cassette subfamily B member                      |
    | TAP2     | transporter 2, ATP binding cassette subfamily B member                      |
    | CD8A     | CD8 subunit alpha                                                           |
    | CD8B     | CD8 subunit beta                                                            |
    | CD8B2    | CD8B family member 2                                                        |
    | KIR3DL2  | killer cell immunoglobulin like receptor, three Ig domains and long cytoplasmic tail 2 |
    | KIR3DL1  | killer cell immunoglobulin like receptor, three Ig domains and long cytoplasmic tail 1 |
    | KIR3DL3  | killer cell immunoglobulin like receptor, three Ig domains and long cytoplasmic tail 3 |
    | KIR3DS1  | killer cell immunoglobulin like receptor, three Ig domains and short cytoplasmic tail 1 |
    | KIR2DL2  | killer cell immunoglobulin like receptor, two Ig domains and long cytoplasmic tail 2 |
    | KIR2DL1  | killer cell immunoglobulin like receptor, two Ig domains and long cytoplasmic tail 1 |
    | KIR2DL3  | killer cell immunoglobulin like receptor, two Ig domains and long cytoplasmic tail 3 |
    | KIR2DL4  | killer cell immunoglobulin like receptor, two Ig domains and long cytoplasmic tail 4 |
    | KIR2DL5A | killer cell immunoglobulin like receptor, two Ig domains and long cytoplasmic tail 5A |
    | KLRC1    | killer cell lectin like receptor C1                                         |
    | KLRC2    | killer cell lectin like receptor C2                                         |
    | KLRC3    | killer cell lectin like receptor C3                                         |
    | KLRC4    | killer cell lectin like receptor C4                                         |
    | KLRD1    | killer cell lectin like receptor D1                                         |
    | KIR2DS1  | killer cell immunoglobulin like receptor, two Ig domains and short cytoplasmic tail 1 |
    | KIR2DS3  | killer cell immunoglobulin like receptor, two Ig domains and short cytoplasmic tail 3 |
    | KIR2DS4  | killer cell immunoglobulin like receptor, two Ig domains and short cytoplasmic tail 4 |
    | KIR2DS5  | killer cell immunoglobulin like receptor, two Ig domains and short cytoplasmic tail 5 |
    | KIR2DS2  | killer cell immunoglobulin like receptor, two Ig domains and short cytoplasmic tail 2 |
    | IFI30    | IFI30 lysosomal thiol reductase                                              |
    | LGMN     | legumain                                                                    |
    | CTSB     | cathepsin B                                                                 |
    | HLA-DMA  | major histocompatibility complex, class II, DM alpha                        |
    | HLA-DMB  | major histocompatibility complex, class II, DM beta                         |
    | HLA-DOA  | major histocompatibility complex, class II, DO alpha                        |
    | HLA-DOB  | major histocompatibility complex, class II, DO beta                         |
    | HLA-DPA1 | major histocompatibility complex, class II, DP alpha 1                      |
    | HLA-DPB1 | major histocompatibility complex, class II, DP beta 1                       |
    | HLA-DQA1 | major histocompatibility complex, class II, DQ alpha 1                      |
    | HLA-DQA2 | major histocompatibility complex, class II, DQ alpha 2                      |
    | HLA-DQB1 | major histocompatibility complex, class II, DQ beta 1                       |
    | HLA-DRA  | major histocompatibility complex, class II, DR alpha                        |
    | HLA-DRB1 | major histocompatibility complex, class II, DR beta 1                       |
    | HLA-DRB3 | major histocompatibility complex, class II, DR beta 3                       |
    | HLA-DRB4 | major histocompatibility complex, class II, DR beta 4                       |
    | HLA-DRB5 | major histocompatibility complex, class II, DR beta 5                       |
    | CD74     | CD74 molecule                                                               |
    | CTSL     | cathepsin L                                                                 |
    | CTSS     | cathepsin S                                                                 |
    | CD4      | CD4 molecule                                                                |
    | CIITA    | class II major histocompatibility complex transactivator                    |
    | RFX5     | regulatory factor X5                                                        |
    | RFXANK   | regulatory factor X associated ankyrin containing protein                   |
    | RFXAP    | regulatory factor X associated protein                                      |
    | CREB1    | cAMP responsive element binding protein 1                                   |
    | NFYA     | nuclear transcription factor Y subunit alpha                                |
    | NFYB     | nuclear transcription factor Y subunit beta                                 |
    | NFYC     | nuclear transcription factor Y subunit gamma                                | 

## Generate a custom interval set

To build a custom set of intervals, you'll need `bedtools` installed. If you built a conda environment for the setup using the `setup_env.yml` file, `bedtools` will be present.

Download from the GIAB ftp the intervals to use

```sh
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/LowComplexity/GRCh38_AllTandemRepeats_51to200bp_slop5.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/LowComplexity/GRCh38_AllTandemRepeats_ge101bp_slop5.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/OtherDifficult/GRCh38_KIR.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/OtherDifficult/GRCh38_MHC.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/OtherDifficult/GRCh38_VDJ.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_quadTR_50to149_slop5.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_quadTR_ge150_slop5.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_triTR_50to149_slop5.bed.gz &\
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.4/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_homopolymer_ge12_slop5.bed.gz
```

These intervals have a slop of 5bp, in both start and end. So we're gonna remove it before intersecting intervals with bedtools with the following script, that we save inside a file called `drop_slop.py` 

```python
import sys
import gzip
import os

def main(iput: str):
    with gzip.open(iput, 'rt') as f:
        outfile = os.path.basename(iput).replace('_slop5', '')
        with gzip.open(outfile, 'wt') as out:
            for line in f:
                chrom, start, end = line.rstrip().split('\t')
                if int(end) >= int(start) and int(end) - int(start) > 1:
                    chrom = chrom.split('chr')[1]
                    start = int(start) + 4
                    end = int(end) - 4
                    if end - start > 5:
                        out.write(f'{chrom}\t{start}\t{end}\n')
                else:
                    continue

if __name__ == '__main__':
    main(sys.argv[1])
```

Then we'll edit all the intervals with

```sh
for bed_file in $PWD/*slop5*.bed.gz ; do python3 drop_slop.py $bed_file ${bed_file/_slop5/_no_slop} ; done
```

Now we could intersect all the intervals with the custom interval file. Assuming it is called `calling_intervals.BED.gz`, the command will be

```sh
files=($(find $PWD -maxdepth 1 -type f -name "GRCh38*.bed.gz" | grep -v "_slop5")) &\
bedtools multiinter -i calling_intervals.BED.gz "${files[@]}" | awk -F'\t' '{if ($5 == 1) print $0}' - | sort -k1,1 -k2,2n | bgzip -c > ENEO_callset.BED.gz
```

!!! tip
    Focusing on exonic portions of target genes is less likely to produce systematic errors in variant calling from RNA-seq, due to the splicing mechanism. A nice trick to get these regions could be using just awk and the gtf:
    ```
    zgrep "protein_coding" gencode.gtf.gz | awk -F'\t' '{ if ($3 == "exon") print $1, $2, $3}' OFS='\t' - | bgzip -c > calling_intervals.BED.gz   
    ```

