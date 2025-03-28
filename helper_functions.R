#####----------------------------------------------------------------------#####
##### generate DESeq2 dataset object
#####----------------------------------------------------------------------#####
generate_DESeq2_dataset <- function(path.to.salmon.quants, 
                                    path.to.tx2gene,
                                    meta.data, 
                                    file.ext = "quant.sf"){
  # fetch all quant.sf files
  files <- file.path(path.to.salmon.quants, meta.data$sample, file.ext)
  names(files) <- meta.data$sample
  
  # read in the tx2gene list
  tx2gene <- read_tsv(path.to.tx2gene, 
                      col_names = c("transcript_id", "gene_id", "gene_name"), 
                      show_col_types = FALSE)
  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
  output <- DESeqDataSetFromTximport( txi,
                                      colData = meta.data,
                                      design = ~ condition) # make sure that group exists in meta.data
  return(output)
}

#####----------------------------------------------------------------------#####
##### Generate interactive table in html file
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####----------------------------------------------------------------------#####
##### RUN DESEQ2
#####----------------------------------------------------------------------#####
run_DESeq2_and_preprocess <- function(  deseq.dataset, 
                                        tx2gene, 
                                        thresh.pval = 0.05,
                                        geneid.to.remove = NULL,
                                        controlGenes = NULL){
  
  # run DESeq2
  if (is.null(geneid.to.remove) == FALSE){
    deseq.dataset <- deseq.dataset[rownames(deseq.dataset) %in% geneid.to.remove == FALSE]
  }
  
  if (is.null(controlGenes) == FALSE){
    print("ERCC spike in, specify in controlGenes when estimating size factors ...")
    deseq.dataset <- estimateSizeFactors(deseq.dataset, controlGenes = controlGenes)
  }
  deseq.object <- DESeq(deseq.dataset)
  
  # normalized read counts
  norm.counts <- counts(deseq.object, normalized=TRUE)
  
  # preprocessing
  ## extract the Differential expression analysis results from DESeq2 object
  ## this object contains: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  deseq.result <- results(deseq.object)
  
  # generate a dataframe containing gene.id and gene.name
  genedf <- data.frame(gene_id = tx2gene$gene_id, 
                       gene_name = tx2gene$gene_name) 
  
  ## remove duplicated gene.name / gene.id
  genedf <- genedf[!duplicated(genedf), ]
  
  ## merge the normalized count matrix and the gene information altogether
  norm.counts <- merge(genedf, norm.counts, by.x= "gene_id", by.y= "row.names")
  
  ## merge the normalized count matrix (with gene information) with DESeq2 results
  resultdf <- merge(genedf, as.data.frame(deseq.result), by.x= "gene_id", by.y="row.names", all.x=F, all.y=T)
  
  resultdf <- merge(resultdf, norm.counts, by=c("gene_id", "gene_name"), all.x=T, all.y=F)
  
  ## remove rows with NA values
  resultdf <- resultdf[complete.cases(resultdf), ]
  
  ## marking significant and non-significant genes by the DESeq2 p-value
  resultdf$sig <- "Non-sig."
  
  resultdf$sig[resultdf$padj < thresh.pval] <- "Sig. genes"
  
  ## Marking ERCC genes
  sel_ERCC <- str_detect(resultdf$gene_id, "^ERCC-*gene")
  
  resultdf$sig[sel_ERCC] <- "Spike in"
  
  ## remove ERCC spike-in genes
  resultdf <- resultdf[!sel_ERCC,]
  
  ## extract the dataframe with significant genes only
  resultdf.sig <- resultdf[resultdf$padj < thresh.pval,]
  
  ## extract the dataframe with significant genes only
  resultdf.nonsig <- resultdf[resultdf$padj >= thresh.pval,]
  
  # collect all outputs
  output <- list(norm.count = norm.counts,
                 all.resultdf = resultdf,
                 resultdf.sig = resultdf.sig,
                 resultdf.nonsig = resultdf.nonsig,
                 deseq.object = deseq.object)
  return(output)
}
