perform_gsea_analysis <- function(counts, condition_list, ref_group, compare_group) {
  # Load required libraries
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  
  # Ensure counts is a matrix
  counts_matrix <- as.matrix(counts)
  
  # Match coldata to counts
  coldata <- condition_list[colnames(counts_matrix), , drop = FALSE]
  
  # Subset to desired groups
  subset_samples <- rownames(coldata)[coldata$group %in% c(compare_group, ref_group)]
  counts_subset <- counts_matrix[, subset_samples]
  coldata_subset <- coldata[subset_samples, , drop = FALSE]
  coldata_subset$group <- factor(coldata_subset$group, levels = c(compare_group, ref_group))
  
  # DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                                colData = coldata_subset,
                                design = ~ group)
  
  # Filter low count genes
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", compare_group, ref_group))
  res <- res[!is.na(res$stat), ]
  
  # Map to Entrez IDs
  gene_symbols <- rownames(res)
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  res$entrez <- entrez_ids
  res <- res[!is.na(res$entrez), ]
  
  # Create ranked list
  gene_list <- res$stat
  names(gene_list) <- res$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Perform GSEA
  gsea_go <- gseGO(geneList = gene_list,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   keyType = "ENTREZID",
                   nPerm = 1000,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = TRUE)
  
  # Return result dataframe
  result_df <- as.data.frame(gsea_go@result)
  return(result_df)
}
