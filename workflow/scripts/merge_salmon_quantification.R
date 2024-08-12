
library("dplyr")
library("tximport")
library("tximeta")
library("tibble")


files <- file.path(snakemake@input[["quant"]])
patients <- snakemake@params[["patients"]]
coldata <- data.frame(files=files, names=patients)
# # build linked Txome
indexDir <- file.path(snakemake@input[["index"]])
fasta <- file.path(snakemake@input[["cdna_fasta"]])
gtf <- file.path(snakemake@input[["annotation"]])
makeLinkedTxome(indexDir=indexDir,source="Ensembl",
organism="Homo sapiens",release="105",genome="GRCh38",
fasta=fasta,gtf=gtf)
se <- tximeta(coldata, useHub=FALSE)

outdir <- file.path(snakemake@params[["outdir"]])

transcript_level_TPM <- se@assays@data@listData[["abundance"]] %>% as.data.frame() %>% tibble::rownames_to_column(var="transcripts")
write.table(transcript_level_TPM, sep="\t", file=paste0(outdir,"/transcript_expression.tsv"),row.names = F, quote=F)
gse <- summarizeToGene(se, countsFromAbundance="lengthScaledTPM")
gene_level_TPM <- gse@assays@data@listData[["abundance"]] %>% as.data.frame() %>% tibble::rownames_to_column(var="genes")
write.table(gene_level_TPM, sep="\t", file=paste0(outdir,'/gene_expression.tsv'),row.names = F, quote=F)

