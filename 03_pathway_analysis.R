gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(comprehenr)
library(ggrepel)
library(heatmaply)
library(org.Mm.eg.db)

if ("org.Mm.eg.db" %in% installed.packages() == FALSE){
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

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
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

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

# i <- 1

for (i in seq(1, length(all.comparisons)) ){
  comp <- all.comparisons[[i]]
  condition1 <- comp[[1]]
  condition2 <- comp[[2]]
  
  path.to.save.output <- file.path(path.to.03.output, sprintf("pathway_analysis_%s_vs_%s", condition1, condition2))
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.deseq.output <- file.path(path.to.02.output, sprintf("DGEA_%s_vs_%s", condition1, condition2))
  deseq.output <- readRDS(file.path(path.to.deseq.output, "deseq_output.rds"))
  
  resdf <- deseq.output$all.resultdf
  
  ######------------------------------------------------------------------------------#####
  ###### INSTALL the newest clusterProfiler from bioconductor
  ######------------------------------------------------------------------------------#####
  # https://bioconductor.org/packages/release/bioc/src/contrib/clusterProfiler_4.14.6.tar.gz
  if (packageVersion("clusterProfiler") != "4.14.6"){
    print("install clusterProfiler 4.14.6 from Bioconductor")
    install.packages("https://cran.r-project.org/src/contrib/yulab.utils_0.2.0.tar.gz", 
                     type = "source", repos = NULL)
    install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/clusterProfiler_4.14.6.tar.gz", 
                     type = "source", repos = NULL)
  } else {
    print(sprintf("clusterProfiler exists, version %s", packageVersion("clusterProfiler/")))
  }
  
  ######------------------------------------------------------------------------------#####
  
  library(clusterProfiler)
  convertdf <- bitr(resdf$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  resdf <- merge(resdf, convertdf, by.x = "gene_name", by.y = "SYMBOL")
  
  ##### MAIN RUN: PATHWAY ANALYSIS
  
  #####--------------------------------------------------------------------#####
  ##### ORA WITH GO
  #####--------------------------------------------------------------------#####
  ora.GO <- enrichGO(gene = subset(resdf, resdf$padj <= 0.05)$gene_name,
                     OrgDb = org.Mm.eg.db,
                     ont = "ALL",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE,
                     keyType = "SYMBOL",
                     pAdjustMethod = "BH")
  
  if (is.null(ora.GO) == TRUE){
    ora.GOdf <- data.frame(status = c("No results from ORA - GO"))  
  } else {
    ora.GOdf <- as.data.frame(ora.GO)
    ora.GOdf <- ora.GOdf[order(ora.GOdf$p.adjust, decreasing = FALSE),]  
  }
  
  writexl::write_xlsx(ora.GOdf, file.path(path.to.03.output, "ORA_GO_result.xlsx"))
  saveRDS(ora.GO, file.path(path.to.save.output, "ora.GO.rds"))
  
  #####--------------------------------------------------------------------#####
  ##### ORA WITH KEGG
  #####--------------------------------------------------------------------#####
  ora.KEGG <-  enrichKEGG(gene = subset(resdf, resdf$padj <= 0.05)$ENTREZID,
                          organism     = 'mmu',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05)
  if (is.null(ora.KEGG) == TRUE){
    ora.KEGGdf <- data.frame(status = c("No results obtained from KEGG-ORA"))
  } else {
    ora.KEGGdf <- as.data.frame(ora.KEGG)
    ora.KEGGdf <- ora.KEGGdf[order(ora.KEGGdf$p.adjust), ]    
  }
  
  writexl::write_xlsx(ora.KEGGdf, file.path(path.to.03.output, "ORA_KEGG_result.xlsx"))
  saveRDS(ora.KEGG, file.path(path.to.save.output, "ora.KEGG.rds"))
  
  #####--------------------------------------------------------------------#####
  ##### GSEA WITH GO
  #####--------------------------------------------------------------------#####
  tmp.full.list <- resdf %>% arrange(desc(log2FoldChange))
  input.gene.list <- tmp.full.list$log2FoldChange
  names(input.gene.list) <- tmp.full.list$gene_name
  
  GSEA.GO <- gseGO(geneList = input.gene.list,
                   OrgD = org.Mm.eg.db,
                   ont = "ALL",
                   minGSSize = 100,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   keyType = "SYMBOL", seed = TRUE)
  
  GSEA.GOdf <- as.data.frame(GSEA.GO) 
  
  GSEA.GOdf <- GSEA.GOdf %>% rowwise() %>% 
    mutate(abs.NES = abs(NES)) %>%
    rownames_to_column("idx")
  
  GSEA.GOdf <- GSEA.GOdf[order(GSEA.GOdf$NES, decreasing = TRUE), ]
  saveRDS(object = GSEA.GO, file.path(path.to.03.output, "GSEA.GO.rds"))
  writexl::write_xlsx(GSEA.GOdf, file.path(path.to.save.output, "GSEA_GOdf.xlsx"))
  
  #####--------------------------------------------------------------------#####
  ##### GSEA WITH KEGG
  #####--------------------------------------------------------------------#####
  tmp.full.list <- resdf %>% arrange(desc(log2FoldChange))
  input.gene.list <- tmp.full.list$log2FoldChange
  names(input.gene.list) <- convertdf$ENTREZID
  
  GSEA.KEGG <- gseKEGG(geneList = input.gene.list,
                       organism = "mmu",
                       minGSSize = 100,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       verbose = FALSE, seed = TRUE)
  
  GSEA.KEGGdf <- as.data.frame(GSEA.KEGG)
  
  GSEA.KEGGdf <- GSEA.KEGGdf %>% rowwise() %>% 
    mutate(abs.NES = abs(NES)) %>%
    rownames_to_column("idx")
  
  GSEA.KEGGdf <- GSEA.KEGGdf[order(GSEA.KEGGdf$NES, decreasing = TRUE), ]
  
  saveRDS(object = GSEA.KEGG, file.path(path.to.save.output, "GSEA.KEGG.rds"))
  writexl::write_xlsx(GSEA.KEGGdf, file.path(path.to.save.output, "GSEA_KEGGdf.xlsx")) 
}

