# Setup resources

ENEO heavily depends on public genetic databases for germline probability estimation. It's required to download data from multiple resources and convert them in order to have concordant annotations (i.e. Ensembl).

## Automatic setup
On working

![](https://imgs.xkcd.com/comics/automation.png)

## Manual setup

To prepare input files, it's required to have an environment with these tools installed 

- bedtools
- bcftools
- tabix
- gatk4

The easiest way is to create a single conda environment as 

```
conda create --name eneo_setup -c bioconda -c conda-forge bedtools bcftools tabix gatk4
```


### Genome and Transcriptome files

Given the use of Ensembl annotation, files could be downloaded from the Ensembl ftp site

**Genome**:
```
wget https://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.dict
```
**Transcriptome**
```
wget https://ftp.ensembl.org/pub/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

### Genetic population resources

ENEO uses multiple databases to infer germline likelihood of candidate variants. Most of them uses different chromosome naming, so you need to convert them. A conversion table is available at

#### dbSNP

```bash
wget -c https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz
wget -c https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz.tbi
```

#### 1000G

```bash
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```
#### known indels

```bash
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
```

### gnomAD

```bash
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
```





