gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(comprehenr)
library(ggrepel)
library(heatmaply)

if ("tximport" %in% installed.packages() == FALSE){
  BiocManager::install("tximport", update = FALSE)
}

library(tximport)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/241203_Eswaran_RNAseq"
source(file.path(path.to.main.src, "helper_functions.R"))
source(file.path(path.to.main.src, "list_of_comparisons.R"))

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "241203_Eswaran_Schippers"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.input <- "/media/hieunguyen/HD01/storage/241203_Eswaran_Schippers_KJM_3mRNAseq_2_Processed_data"
path.to.star.salmon <- file.path(path.to.main.input, "nfcore_3mRNAseq/results/star_salmon")
path.to.tx2gene <- file.path(path.to.star.salmon, "tx2gene.tsv")
tx2gene <- read_tsv(path.to.tx2gene, 
                    col_names = c("transcript_id", "gene_id", "gene_name"), 
                    show_col_types = FALSE)

all.sample.names <- Sys.glob(file.path(path.to.star.salmon, "*.markdup.sorted.bam")) %>% basename()
all.sample.names <- to_vec(for (item in all.sample.names) str_replace(item, ".markdup.sorted.bam", "")) %>% unique()

meta.data <- data.frame(
  sample = all.sample.names,
  condition = to_vec( for (item in all.sample.names){
    paste0(str_split(item, "_")[[1]][1:length(str_split(item, "_")[[1]]) - 1], collapse = "_")
  } )
)
meta.data$condition <- factor(meta.data$condition, levels = unique(meta.data$condition))

for (i in seq(1, length(all.comparisons)) ){
  comp <- all.comparisons[[i]]
  condition1 <- comp[[1]]
  condition2 <- comp[[2]]
  
  print(sprintf("Working on differential gene expression analysis between %s vs %s", condition1, condition2))
  
  path.to.save.output <- file.path(path.to.02.output, sprintf("DGEA_%s_vs_%s", condition1, condition2))
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  
  filtered.metadata <- subset(meta.data, meta.data$condition %in% c(condition1, condition2))
  filtered.metadata$condition <- factor(filtered.metadata$condition, levels = c(condition1, condition2))
  writexl::write_xlsx(filtered.metadata, file.path(path.to.save.output, "metadata.xlsx"))
  
  condition1.samples <- subset(filtered.metadata, filtered.metadata$condition == condition1)$sample
  condition2.samples <- subset(filtered.metadata, filtered.metadata$condition == condition2)$sample
  
  deseq.dataset <-generate_DESeq2_dataset(path.to.star.salmon, 
                                          path.to.tx2gene = path.to.tx2gene,
                                          meta.data = filtered.metadata, 
                                          file.ext = "quant.sf")
  
  ##### get ERCC spike in genes
  all.genes <- row.names(deseq.dataset)
  ERCC.genes <- to_vec( for (item in all.genes) {
    if (grepl("ERCC", item) == TRUE){
      item
    }
  })
  
  ##### run deseq main analysis
  deseq.output <- run_DESeq2_and_preprocess(deseq.dataset = deseq.dataset, 
                                            tx2gene = tx2gene, 
                                            thresh.pval = 0.05,
                                            controlGenes = which(row.names(deseq.dataset) %in% ERCC.genes))
  ##### save intermediate data file deseq.output
  saveRDS(deseq.output, file.path(path.to.save.output, "deseq_output.rds"))
  
  ##### export table of significantly expressed genes
  sigdf <- deseq.output$resultdf.sig %>%
    rowwise() %>% 
    mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
    arrange(desc(abs_log2FoldChange))
  
  sigdf[[sprintf("%s_baseMean", condition1)]] <- rowMeans(sigdf[, c(condition1.samples)])
  sigdf[[sprintf("%s_baseMean", condition2)]] <- rowMeans(sigdf[, c(condition2.samples)])
  writexl::write_xlsx(sigdf, file.path(path.to.save.output, "sig_genes.xlsx"))
  
  ##### export table of nonsignificantly expressed genes
  nonsigdf <- deseq.output$resultdf.nonsig
  nonsigdf[[sprintf("%s_baseMean", condition1)]] <- rowMeans(nonsigdf[, c(condition1.samples)])
  nonsigdf[[sprintf("%s_baseMean", condition2)]] <- rowMeans(nonsigdf[, c(condition2.samples)])
  writexl::write_xlsx(nonsigdf, file.path(path.to.save.output, "nonsig_genes.xlsx"))
  
  ##### PCA plot
  input.df <- deseq.output$norm.count[filtered.metadata$sample]
  pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)
  
  pcadf <- data.frame(pca.object$x)
  row.names(pcadf) <- filtered.metadata$sample
  
  pcadf <- merge(pcadf, filtered.metadata, by.x = "row.names", by.y = "sample")
  
  pca.plot <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
    geom_point(size = 4) +
    ggtitle(sprintf("PCA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    geom_label_repel() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
  
  ggsave(plot = pca.plot, filename = "PCA.svg", 
         path = path.to.save.output, 
         dpi = 300, 
         width = 10, 
         height = 10)
  
  ##### Volcano plot
  cutoff.adjp <- 0.05
  cutoff.logFC <- 1
  
  input.df <- deseq.output$all.resultdf
  input.df <- input.df %>% mutate(abs.log2FoldChange = abs(log2FoldChange))
  input.df <- input.df %>% rowwise() %>%
    mutate(show.gene.name = ifelse(padj < cutoff.adjp, gene_name, NA))
  
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=log2FoldChange, y=-log10(padj), col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    geom_text_repel() +
    ggtitle(sprintf("Volcano plot, fold-change = %s / %s (right/left)", condition2, condition1)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    xlim(c(-max(input.df$abs.log2FoldChange), max(input.df$abs.log2FoldChange)))
  ggsave(plot = volcano.plot, filename = "volcano_plot.svg", path = path.to.save.output, dpi = 300, width = 14, height = 10)
  
  ##### MA plot
  ma.plot <- ggplot(data=input.df, 
                    aes(x=log2(baseMean), y=log2FoldChange, col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_text_repel() +
    ggtitle(sprintf("MA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12))
  ggsave(plot = ma.plot, filename = "MA_plot.svg", path = path.to.save.output, dpi = 300, width = 14, height = 10)
  
  ##### Heatmap
  sig.genes.with.highlogFC <- subset(input.df, (input.df$sig == "Sig. genes") & (input.df$abs.log2FoldChange > cutoff.logFC))
  nonsig.genes.with.highlogFC <- subset(input.df, (input.df$sig != "Sig. genes") & (input.df$abs.log2FoldChange > cutoff.logFC))
  input.to.heatmap <- subset(sig.genes.with.highlogFC, select = c("gene_id", "gene_name", filtered.metadata$sample))
  
  if (nrow(input.to.heatmap) == 1){
    p <- ggplot() + ggtitle("No or only 1 significantly differently expressed genes, cannot show heatmap")
    ggsave(plot = p, filename = "no_data_to_show_heatmap.svg", path = file.path(path.to.save.output), device = "svg", width = 14, height = 10, dpi = 300)
  } else {
    if (nrow(input.to.heatmap) > 0){
      heatmap.values <- log10(input.to.heatmap[,3:(dim(input.to.heatmap)[2])] + 1)
      selected.genes.heatmap.plot <- heatmaply(heatmap.values, 
                                               main= sprintf("Heatmap, Sample: %s vs. %s", condition1, condition2),
                                               method = "plotly",labRow=input.to.heatmap$gene_name,
                                               xlab = "Samples", ylab = "Genes", width = 800, height = 600,
                                               showticklabels = c(TRUE, FALSE), show_dendrogram = c(FALSE, TRUE),
                                               key.title = "log10 scale colormap",
                                               label_names = c("Gene", "Sample", "Expression"),
                                               k_col = 2, file = file.path(path.to.save.output, "heatmap.html"))
    }
  }
  print(sprintf("Finish generating results for DGEA %s vs %s", condition1, condition2))
}
