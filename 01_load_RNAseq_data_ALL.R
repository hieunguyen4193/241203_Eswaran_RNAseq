gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(comprehenr)
library(ggrepel)
if ("tximport" %in% installed.packages() == FALSE){
  BiocManager::install("tximport", update = FALSE)
}

library(tximport)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/241203_Eswaran_RNAseq"
source(file.path(path.to.main.src, "helper_functions.R"))

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "241203_Eswaran_Schippers"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.input <- "/media/hieunguyen/HD01/storage/241203_Eswaran_Schippers_KJM_3mRNAseq_2_Processed_data"
path.to.star.salmon <- file.path(path.to.main.input, "nfcore_3mRNAseq/results/star_salmon")
path.to.tx2gene <- file.path(path.to.star.salmon, "tx2gene.tsv")

all.sample.names <- Sys.glob(file.path(path.to.star.salmon, "*.markdup.sorted.bam")) %>% basename()
all.sample.names <- to_vec(for (item in all.sample.names) str_replace(item, ".markdup.sorted.bam", "")) %>% unique()

meta.data <- data.frame(
  sample = all.sample.names,
  condition = to_vec( for (item in all.sample.names){
    paste0(str_split(item, "_")[[1]][1:length(str_split(item, "_")[[1]]) - 1], collapse = "_")
  } )
)

meta.data$condition <- factor(meta.data$condition, levels = unique(meta.data$condition))

deseq.dataset <- generate_DESeq2_dataset(path.to.star.salmon, 
                                         path.to.tx2gene = path.to.tx2gene,
                                         meta.data = meta.data, 
                                         file.ext = "quant.sf")
tx2gene <- read_tsv(path.to.tx2gene, 
                    col_names = c("transcript_id", "gene_id", "gene_name"), 
                    show_col_types = FALSE)

##### apply variance stabilizing transformations (VST) 
vsd <- vst(deseq.dataset, blind=FALSE)
exprs.matrix <- assay(vsd)

all.genes <- row.names(exprs.matrix)
ERCC.genes <- to_vec( for (item in all.genes) {
  if (grepl("ERCC", item) == TRUE){
    item
  }
})

deseq.dataset.noERCC <- deseq.dataset[!rownames(deseq.dataset) %in% ERCC.genes,]
vsd.no.ERCC <- vst(deseq.dataset.noERCC, blind=FALSE)
exprs.matrix.no.ERCC <- assay(vsd.no.ERCC)

##### generate PCA for all data
pca.object <- prcomp(t(exprs.matrix), rank. = 2, scale. = FALSE)
pcadf <- data.frame(pca.object$x)
pcadf <- merge(pcadf, meta.data, by.x = "row.names", by.y = "sample")

# PCA of all samples
pca.plot.all.samples <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
  geom_point(size = 4) +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
  geom_text_repel(nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 20)+ 
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  xlim(c(-150,150)) + ylim(c(-150, 150)) 

# assess replicates
pca.plot.replicates <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
  geom_point(size = 4) +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
  geom_text_repel(nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 20)+ 
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  xlim(c(-150,150)) + ylim(c(-150, 150)) + 
  facet_wrap(~condition)

##### generate PCA for all data
pca.object <- prcomp(t(exprs.matrix), rank. = 2, scale. = FALSE)
pcadf <- data.frame(pca.object$x)
pcadf <- merge(pcadf, meta.data, by.x = "row.names", by.y = "sample")

# PCA of all samples
pca.plot.all.samples.noERCC <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
  geom_point(size = 4) +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
  geom_text_repel(nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 20)+ 
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  xlim(c(-150,150)) + ylim(c(-150, 150)) 

# assess replicates
pca.plot.replicates.noERCC <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
  geom_point(size = 4) +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
  geom_text_repel(nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 20)+ 
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  xlim(c(-150,150)) + ylim(c(-150, 150)) + 
  facet_wrap(~condition)




