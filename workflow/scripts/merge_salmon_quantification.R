library("dplyr")
library("tximport")
library("rtracklayer")
library("tibble")

files <- file.path(snakemake@input[["quant"]])
patients <- snakemake@params[["patients"]]
gtf_file <- file.path(snakemake@input[["annotation"]])

# Label each quantification file with its patient so the output columns carry
# the sample identity instead of anonymous defaults.
names(files) <- patients

gtf <- rtracklayer::import(gtf_file)

tx2gene <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_id")])
tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]
# The GTF carries one row per feature (exon/CDS/...), so each transcript->gene
# mapping is repeated many times; keep the unique pairs only.
tx2gene <- unique(tx2gene)

# Transcript-level import: txOut = TRUE keeps the per-transcript estimates
# instead of summarising them to gene level.
txi_tx <- tximport(files,
                   type = "salmon",
                   txOut = TRUE,
                   ignoreTxVersion = TRUE,
                   countsFromAbundance = "lengthScaledTPM")

transcript_level_TPM <- txi_tx$abundance %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "transcripts")

write.table(transcript_level_TPM,
            file = snakemake@output[["transcript"]],
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Summarise the transcript-level estimates to gene level.
txi_gene <- summarizeToGene(txi_tx,
                            tx2gene = tx2gene,
                            ignoreTxVersion = TRUE,
                            countsFromAbundance = "lengthScaledTPM")

gene_level_TPM <- txi_gene$abundance %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "genes")

write.table(gene_level_TPM,
            file = snakemake@output[["gene"]],
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
