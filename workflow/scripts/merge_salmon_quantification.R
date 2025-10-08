library("dplyr")
library("tximport")
library("rtracklayer") 
library("tibble")

files <- file.path(snakemake@input[["quant"]])
patients <- snakemake@params[["patients"]]
gtf_file <- file.path(snakemake@input[["annotation"]])
coldata <- data.frame(files = files, names = patients, stringsAsFactors = FALSE)

gtf <- rtracklayer::import(gtf_file)

tx2gene <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_id")])

tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE,
                countsFromAbundance = "lengthScaledTPM")
transcript_level_TPM <- txi$abundance %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "transcripts")

write.table(transcript_level_TPM, 
            file = snakemake@output[["transcript"]],
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

gene_level_TPM <- txi$abundance %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "genes")

write.table(gene_level_TPM, 
            file = snakemake@output[["gene"]],
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)