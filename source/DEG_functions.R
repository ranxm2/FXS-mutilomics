library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(pheatmap)
library(tibble)
library(rtracklayer)
library(tidyr)
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(gprofiler2)
library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)
library(edgeR)
library(parallel)
library(BiocParallel)
library(GSVA)
library(decoupleR)
library(ggrepel)
library(patchwork)
library(ggsignif) 
library(extrafont)
library(GSEABase)
library(ggsignif)

Result_folder_structure <- function(result_folder) {
  dir.create(result_folder, showWarnings = FALSE)
  dir.create(paste0(result_folder, "/01-Sample_info"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/02-DEG"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/03-Enrichment"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/04-GSVA"), showWarnings = FALSE)
  # dir.create(paste0(result_folder, "/05-22q_Gene"), showWarnings = FALSE)
}

# Test for the function
# result_folder = "results/02-DEG-Vehicle"
# Result_folder_structure(result_folder)



DEAnalysis <- function(counts =NULL, 
                       reference_group = NULL, 
                       compare_group = NULL, 
                       condition_list = NULL, 
                       target_gene =NULL, 
                       result_folder = NULL) {
  # Filter the count data 
  sample_info <- condition_list %>%
    mutate(sample = rownames(.)) %>%
    filter(group %in% c(reference_group, compare_group))
  count_input <- counts[, sample_info$sample]
  
  # filter gene withall 0
  count_input <- count_input[rowSums(count_input) > 0, ]
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_input,
                                colData = sample_info,
                                design = ~ group)
  dds$group <- relevel(dds$group, ref = reference_group)
  

  ##############################
  ### Run DESeq analysis    ####
  ###------------------------###
  
  dds <- DESeq(dds)
  return(dds)
}



plot_sample_heatmap <- function(dds_obj,figure_folder, file_name) {
   # Now apply variance stabilizing transformation
   vsd.obj <- varianceStabilizingTransformation(dds_obj, blind = TRUE)
   sampleDists <- dist(t(assay(vsd.obj)))
   sampleDistMatrix <- as.matrix( sampleDists )
   rownames(sampleDistMatrix) <- paste( vsd.obj$group )
   colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
   p <- pheatmap::pheatmap(sampleDistMatrix,
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           col = colors) 
   ggsave(file.path(figure_folder, paste0(file_name, ".png")), 
          p, width = 8, height = 6, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), 
         p, width = 8, height = 6, units = "in", dpi = 300)
   print("Sample distance heatmap is done")
}

# # Test for the function
# plot_sample_heatmap(dds_obj, figure_folder = 
#                       file.path(result_folder,"01-Sample_info"), 
#                     file_name = "01_sample_distance_heatmap")




plot_sample_PCA_plot <- function(dds_obj, figure_folder, file_name, 
                               reference_group, compare_group){

  vsd.obj <- varianceStabilizingTransformation(dds_obj, blind = TRUE)
  pcaData <- plotPCA(vsd.obj,  intgroup = c("group"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Define unique groups and their corresponding colors
  unique_groups <- c(reference_group, compare_group)
  color_palette <- c("#10d5da", "#fe867f")  # Extend this list if you have more groups
  assigned_colors <- setNames(color_palette[seq_along(unique_groups)], unique_groups)
  
  p <-ggplot(pcaData, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
    ) +
    stat_ellipse(level = 0.95)+
    theme_bw() +
    # theme_classic()+
    theme(text = element_text(family = "Arial", colour = "black")) +
    scale_color_manual(values = assigned_colors) +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), 
         p,width = 6, height = 4, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), 
         p,width = 6, height = 4, units = "in")
  print("PCA plot is done")
  print(p)
}

# # test for the function
# plot_sample_PCA_plot(dds_obj,figure_folder = file.path(result_folder,"01-Sample_info"),
#                    file_name = "02_sample_PCA_plot", 
#                    reference_group, compare_group)


plot_volcano_plot <- function(result_df, figure_folder, file_name, 
                              thread = 1, dot_size =2, label_gene = NULL,
                              fig.width = 8, fig.height = 8){
  
  result_df <- result_df  %>% filter(!is.na(padj))
  max_p <- max(-log10(result_df$padj[is.finite(-log10(result_df$padj))]), na.rm = TRUE)

  max_fc <- max(abs(result_df$log2FoldChange), na.rm = TRUE, infinite = FALSE)
  
  
  # Define the factor levels so that "NO" is last
  result_df <- result_df %>%
    mutate(
      diffexpressed = if_else(log2FoldChange > thread & padj < 0.05, "UP",
                              if_else(log2FoldChange < -thread & padj < 0.05, "DOWN", "NO")),
      diffexpressed = factor(diffexpressed, levels = c("NO", "UP", "DOWN" )) # Custom order
    ) %>%
    arrange(diffexpressed)  # Arrange based on the custom factor order
  
  
  label_gene <- intersect(label_gene, result_df$GeneName)
  label_df  <- result_df %>% 
                filter(GeneName %in% label_gene & diffexpressed != "NO")
  
  num_up   <- sum(result_df$diffexpressed == "UP")
  num_down <- sum(result_df$diffexpressed == "DOWN")
  
  p <-ggplot(  result_df , aes(log2FoldChange, -log10(padj), color = diffexpressed)) +
    geom_point(show.legend = F,size =dot_size) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
    scale_x_continuous(limits = c(- max_fc , max_fc )) +
    scale_y_continuous(limits = c(0,max_p  )) +
    # geom_label_repel(max.overlaps = 20, show.legend = F, nudge_y = 2, nudge_x = 2) +
    geom_vline(xintercept = c(-thread, thread)) +
    geom_hline(yintercept = -log10(0.05)) +
    theme_classic() +
    theme(legend.position = "right") +
    annotate("text", x = 0.4 * max_fc + thread, y = 0.9 * max_p, 
             label = paste("UP:", num_up), color = "red",size=6 ) +
    annotate("text", x = -0.4 * max_fc - thread, y = 0.9 * max_p, 
             label = paste("DOWN:", num_down), color = "blue",size=6) +
    
    # Auto-annotations with dynamic colors
    geom_label_repel(data = label_df, 
                     aes(label = GeneName, color = diffexpressed), 
                     size = 5, 
                     max.overlaps =50,
                     min.segment.length = 0,
                     force = 2,
                     show.legend = FALSE, 
                     box.padding = 0.5, 
                     point.padding = 0.3) +
    
    theme( text = element_text(family = "Arial", colour = "black"),
           axis.title = element_text(color = "black", family = "Arial", size = 16),  # Increase axis title size
           axis.text = element_text(color = "black", family = "Arial", size = 14),   # Increase axis text size
           legend.position = "right"
    )
  
  ggsave(file.path(figure_folder, paste0(file_name, ".png")),
         p, width = fig.width, height = fig.height, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")),
         p, width = fig.width, height = fig.height, units = "in")
  print(sprintf("Volcano plot for %s", file_name))
  print(p)
}


# # test for the function
# plot_volcano_plot(result_df, 
#                   figure_folder = file.path(result_folder,"02-DEG"),
#                   file_name = "02_volcano_plot_log2fc_1",
#                   thread = 1 ,label_gene = NULL)
# 
# plot_volcano_plot(result_df, 
#                   figure_folder = file.path(result_folder,"02-DEG"),
#                   file_name = "03_volcano_plot_log2fc_1.5",
#                   thread = 1.5 ,label_gene = NULL)


plot_gene_heatmap <-  function(vsd_obj, gene_list, figure_folder, file_name, 
                               reference_group, compare_group,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               scale = "none",show_rownames =FALSE){
  
  # make the color scale
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  
  # color map
  unique_groups <- c(reference_group, compare_group)
  color_palette <- c("#10d5da","#fe867f")  # Extend this list if you have more groups
  assigned_colors <- setNames(color_palette[seq_along(unique_groups)], unique_groups)
  annotation_colors <- list(group = assigned_colors)
  annotation_col <- as.data.frame(colData(vsd_obj)[, "group", drop = FALSE])

  stabilized_counts <- assay(vsd_obj)
  row_variances <- rowVars(stabilized_counts)
  
  stabilized_counts <- stabilized_counts - rowMeans(stabilized_counts, na.rm=T)
  stabilized_counts <- stabilized_counts[gene_list, rownames(annotation_col)]
  
  print(sprintf("Heatmap for %s ", file_name))
  
  p <-pheatmap::pheatmap(stabilized_counts,
                           color = mr,
                           scale = scale,
                           annotation_col = annotation_col,
                           annotation_colors = annotation_colors,
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols ,
                           fontsize_col = 10,
                           fontsize_row = 10,
                           treeheight_row = 0,
                           show_rownames = show_rownames,
                           border_color = NA)
    
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), 
         p, width = 8, height = 8, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")),
         p, width = 8, height = 8, units = "in")
         
}

# # test for the function
# plot_gene_heatmap(vsd_obj, gene_list = DEG_gene_1, 
#                   figure_folder = file.path(result_folder,"02-DEG"),
#                   file_name = "02_heatmap_log2fc_1_row",
#                   reference_group, compare_group,
#                   cluster_rows = TRUE, cluster_cols = FALSE, 
#                   scale = "row")




map_genes <- function(result, mapping, gene_name_map) {
  # Ensure the result data frame has columns to store the mapped values
  result$intersection_gene <- NA
  result$intersection_gene_id <- NA
  
  for (i in 1:nrow(result)) {
    # Step 1: Extract and split gene IDs from the intersection column
    gene_id <- result$intersection[i]
    gene_id <- strsplit(gene_id, ",")[[1]]
    
    # Step 2: Map gene IDs to Ensembl IDs using the provided mapping
    ensembl_ids <- sapply(gene_id, function(id) mapping[[id]])
    
    # Step 3: Convert the Ensembl IDs to a single string
    ensembl_ids_string <- paste(ensembl_ids, collapse = ",")
    
    # Step 4: Map Ensembl IDs to gene names using the gene_name_map
    gene_names <- sapply(ensembl_ids, function(id) gene_name_map[id])
    gene_names_string <- paste(gene_names, collapse = ",")
    
    # Step 5: Store the results in the result data frame
    result$intersection_gene[i] <- gene_names_string
    result$intersection_gene_id[i] <- ensembl_ids_string
  }
  
  return(result)
}



Enrichment_analysis <- function(gene_list, result_folder, 
                                file_name,gene_name_mapping, flag = "Up"){
  enrichment <- gost(query = gene_list,organism = "hsapiens", correction_method = "fdr", evcodes = T)
  result <- subset(enrichment$result, select = -c(parents))
  result <- map_genes(result, enrichment$meta$genes_metadata$query$query_1$mapping, gene_name_mapping)
  result <-  result %>%  
    mutate(term = sprintf("%s - %s",  source, term_name))
  write.csv(result, file.path(result_folder, paste0(file_name, ".csv")), 
            row.names = FALSE)
  
  print(sprintf("Enrichment analysis for %s ", file_name))
  
  result_df_30 <- result %>%
    dplyr::select(term, p_value) %>%
    arrange(p_value) %>%
    slice_head(n = 30) %>%
    mutate(term = factor(term, levels = rev(unique(term))))
  # Select the top 30 smallest p-values
  
  p<-ggplot(result_df_30, aes(x = term, y = -log10(p_value), fill = -log10(p_value))) +
    geom_bar(stat = "identity", show.legend = TRUE) +
    # Apply different color scales based on flag value
    {
      if (flag == "Up") {
        scale_fill_gradient(high = "#a50f15", low = "#fc9272")  # Red scale for UP
      } else if (flag == "Down") {
        scale_fill_gradient(low = "#56B1F7", high = "#132B43")  # Blue scale for DOWN
      } else {
        scale_fill_manual(values = rep("lightblue", nrow(result_df_30)))  # Default to light blue
      }
    } +
    ylab("-log10(FDR)") +
    xlab("") +
    ggtitle("Top 30 Pathways (Sorted by FDR)") +
    coord_flip() +
    theme_classic(base_size = 14) +
    theme(text = element_text(family = "Arial", size = 14, colour = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(family = "Arial", size = 14, face = "bold", hjust = 0.5)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
  
  print(p)
  ggsave(file.path(result_folder,
                   paste0(file_name, "_top30_pathways.png")), p,
         width = 12, height = 8, units = "in", dpi = 300)
  ggsave(file.path(result_folder, 
                   paste0(file_name, "_top30_pathways.pdf")), p,
         width = 12, height = 8, units = "in")
  
  print(sprintf("Enrichment analysis for GOBP %s ", file_name))
  result_GOBP <- result %>% filter(source == "GO:BP")
  write.csv(result_GOBP, file.path(result_folder, 
                                   paste0(file_name, "_GO_BP.csv")), 
            row.names = FALSE)
  

  result_df_30_GOBP <- result_GOBP %>%
    dplyr::select(term, p_value) %>%
    arrange(p_value) %>%
    slice_head(n = 30) %>%
    mutate(term = factor(term, levels = rev(unique(term))))
  # Select the top 30 smallest p-values
  
  p<-ggplot(result_df_30_GOBP, aes(x = term, y = -log10(p_value), fill = -log10(p_value))) +
    geom_bar(stat = "identity", show.legend = TRUE) +
    # Apply different color scales based on flag value
    {
      if (flag == "Up") {
        scale_fill_gradient(high = "#a50f15", low = "#fc9272")  # Red scale for UP
      } else if (flag == "Down") {
        scale_fill_gradient(low = "#56B1F7", high = "#132B43")  # Blue scale for DOWN
      } else {
        scale_fill_manual(values = rep("lightblue", nrow(result_df_30)))  # Default to light blue
      }
    } +
    ylab("-log10(FDR)") +
    xlab("") +
    ggtitle("Top 30 Pathways for GOBP (Sorted by FDR)") +
    coord_flip() +
    theme_classic(base_size = 14) +
    theme(text = element_text(family = "Arial", size = 14, colour = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(family = "Arial", size = 14, face = "bold", hjust = 0.5)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
  
  print(p)
  ggsave(file.path(result_folder, 
                   paste0(file_name, "_top30_pathways_GOBP.png")), p,
         width = 12, height = 8, units = "in", dpi = 300)
  ggsave(file.path(result_folder, 
                   paste0(file_name, "_top30_pathways_GOBP.pdf")), p,
         width = 12, height = 8, units = "in")
  
  
}




plot_gsva_boxplot<-  function(gsva_matrix, condition_list_label, pathway_name,
                                figure_folder, file_name, 
                                fig.height = 6, fig.width = 4,
                                reference_group, compare_group){
  
  # Sample list sub
  sample_info <- condition_list_label %>%
    mutate(sample = rownames(.)) %>%
    filter(group %in% c(reference_group, compare_group))
  
  # Convert gsva_matrix to a data frame and reshape
  plot_df <- as.data.frame(gsva_matrix) %>%
    rownames_to_column(var = "pathway") %>%
    filter(pathway == pathway_name) %>%  # Select the specific pathway
    pivot_longer(cols = -pathway, names_to = "sample", values_to = "GSVA_score") %>%
    dplyr::select(-pathway)  %>% 
    filter(sample %in% sample_info$sample)  %>%
    left_join(sample_info, by = "sample")   %>% 
    mutate(group = factor(group, levels = c(reference_group, compare_group)))  # Define the factor levels

  
  # Ensure colors are mapped to exact group names
  color_palette <- setNames(c("#10d5da", "#fe867f"), c(reference_group, compare_group))
  
  # Process pathway name: remove first part, capitalize first letter, replace underscores with spaces
  formatted_title <- pathway_name %>%
    str_remove("^[^_]+_") %>%  # Remove everything before the first underscore
    str_to_lower() %>%  # Convert everything to lowercase
    str_replace_all("_", " ") %>%  # Replace underscores with spaces
    str_to_sentence()  # Capitalize first letter
  
  # sub the rna with RNA in formatted_title
  formatted_title <- str_replace_all(formatted_title, "rna", "RNA")
  
  # Create the box plot with scatter overlay
  p<-ggplot(plot_df, aes(x = group, y = GSVA_score)) +
    geom_boxplot(aes(fill = group), alpha = 0.9, outlier.shape = NA, color = "black") +  # Box plot with fill color
    geom_jitter(aes(fill = group), shape = 21, width = 0.2, size = 3, alpha = 0.9, color = "black") +  # Scatter points with fill color
    # stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +  # Mean point
    geom_signif(comparisons = list(c(reference_group, compare_group)), 
                test = "t.test", 
                map_signif_level = TRUE) +  # Add significance annotation
    scale_fill_manual(values = color_palette) +  # Apply custom colors to boxes & scatter dots
    theme_minimal(base_family = "Arial") +  # Use Arial font for all text
    labs(title = formatted_title,  # Use formatted pathway name
         x = "",  # Remove x-axis label
         y = "GSVA Score") +
    theme(text = element_text(size = 14, family = "Arial"),  # Ensure all text uses Arial
          plot.title = element_text(hjust = 0.5, size = 14),  # Center title
          axis.text.x = element_text(size = 16, color = "black"),  # Increase x-axis text size and set color to black
          axis.text.y = element_text(size = 14, color = "black"),  # Set y-axis text color to black
          axis.title.x = element_text(size = 14, color = "black"),  # Set x-axis title color to black
          axis.title.y = element_text(size = 14, color = "black"),  # Set y-axis title color to black
          legend.position = "none")  
  
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), p,
          width = fig.width, height = fig.height, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), p,
         width = fig.width, height = fig.height, units = "in")
  
  return(p)

}


plot_gsva_boxplot_mutil<-  function(gsva_matrix, condition_list_label, pathway_name,
                                      figure_folder, file_name, 
                                      fig.height = 6, fig.width = 8){
  
  # Sample list sub
  sample_info <- condition_list_label %>%
    mutate(sample = rownames(.)) 
  
  # Convert gsva_matrix to a data frame and reshape
  plot_df <- as.data.frame(gsva_matrix) %>%
    rownames_to_column(var = "pathway") %>%
    filter(pathway == pathway_name) %>%  # Select the specific pathway
    pivot_longer(cols = -pathway, names_to = "sample", values_to = "GSVA_score") %>%
    dplyr::select(-pathway)  %>% 
    filter(sample %in% sample_info$sample)  %>%
    left_join(sample_info, by = "sample") 
  
  # Ensure colors are mapped to exact group names
  # color_palette <- setNames(c("#10d5da", "#fe867f"), c(reference_group, compare_group))
  
  # Process pathway name: remove first part, capitalize first letter, replace underscores with spaces
  formatted_title <- pathway_name %>%
    str_remove("^[^_]+_") %>%  # Remove everything before the first underscore
    str_to_lower() %>%  # Convert everything to lowercase
    str_replace_all("_", " ") %>%  # Replace underscores with spaces
    str_to_sentence()  # Capitalize first letter
  
  # sub the rna with RNA in formatted_title
  formatted_title <- str_replace_all(formatted_title, "rna", "RNA")
  
  # Create the box plot with scatter overlay
  p<-ggplot(plot_df, aes(x = group, y = GSVA_score)) +
    geom_boxplot(aes(fill = group), alpha = 0.9, outlier.shape = NA, color = "black") +  # Box plot with fill color
    geom_jitter(aes(fill = group), shape = 21, width = 0.2, size = 3, alpha = 0.9, color = "black") +  # Scatter points with fill color
    # stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +  # Mean point
    # geom_signif(comparisons = list(c(reference_group_1, compare_group_1),
    #                                c(reference_group_2, compare_group_2),
    #                                c(reference_group_3, compare_group_3)),
    #             test = "t.test", 
    #             map_signif_level = TRUE) +  # Add significance annotation
    # scale_fill_manual(values = color_palette) +  # Apply custom colors to boxes & scatter dots
    theme_minimal(base_family = "Arial") +  # Use Arial font for all text
    labs(title = formatted_title,  # Use formatted pathway name
         x = "",  # Remove x-axis label
         y = "GSVA Score") +
    theme(text = element_text(size = 14, family = "Arial"),  # Ensure all text uses Arial
          plot.title = element_text(hjust = 0.5, size = 14),  # Center title
          axis.text.x = element_text(size = 16, color = "black"),  # Increase x-axis text size and set color to black
          axis.text.y = element_text(size = 14, color = "black"),  # Set y-axis text color to black
          axis.title.x = element_text(size = 14, color = "black"),  # Set x-axis title color to black
          axis.title.y = element_text(size = 14, color = "black"),  # Set y-axis title color to black
          legend.position = "none")  
  
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), p,
         width = fig.width, height = fig.height, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), p,
         width = fig.width, height = fig.height, units = "in")
  
  return(p)
  
}




plot_gsva_boxplot_mutil_3<-  function(gsva_matrix, condition_list_label, pathway_name,
                              figure_folder, file_name, 
                              fig.height = 6, fig.width = 8,
                              reference_group_1, compare_group_1,
                              reference_group_2, compare_group_2,
                              reference_group_3, compare_group_3){
  
  # Sample list sub
  sample_info <- condition_list_label %>%
    mutate(sample = rownames(.)) 
  
  # Convert gsva_matrix to a data frame and reshape
  plot_df <- as.data.frame(gsva_matrix) %>%
    rownames_to_column(var = "pathway") %>%
    filter(pathway == pathway_name) %>%  # Select the specific pathway
    pivot_longer(cols = -pathway, names_to = "sample", values_to = "GSVA_score") %>%
    dplyr::select(-pathway)  %>% 
    filter(sample %in% sample_info$sample)  %>%
    left_join(sample_info, by = "sample") 
  
  # Ensure colors are mapped to exact group names
  # color_palette <- setNames(c("#10d5da", "#fe867f"), c(reference_group, compare_group))
  
  # Process pathway name: remove first part, capitalize first letter, replace underscores with spaces
  formatted_title <- pathway_name %>%
    str_remove("^[^_]+_") %>%  # Remove everything before the first underscore
    str_to_lower() %>%  # Convert everything to lowercase
    str_replace_all("_", " ") %>%  # Replace underscores with spaces
    str_to_sentence()  # Capitalize first letter
  
  # sub the rna with RNA in formatted_title
  formatted_title <- str_replace_all(formatted_title, "rna", "RNA")
  
  # Create the box plot with scatter overlay
  p<-ggplot(plot_df, aes(x = group, y = GSVA_score)) +
    geom_boxplot(aes(fill = group), alpha = 0.9, outlier.shape = NA, color = "black") +  # Box plot with fill color
    geom_jitter(aes(fill = group), shape = 21, width = 0.2, size = 3, alpha = 0.9, color = "black") +  # Scatter points with fill color
    # stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +  # Mean point
    geom_signif(comparisons = list(c(reference_group_1, compare_group_1),
                                    c(reference_group_2, compare_group_2),
                                    c(reference_group_3, compare_group_3)),
                test = "t.test", 
                map_signif_level = TRUE) +  # Add significance annotation
    # scale_fill_manual(values = color_palette) +  # Apply custom colors to boxes & scatter dots
    theme_minimal(base_family = "Arial") +  # Use Arial font for all text
    labs(title = formatted_title,  # Use formatted pathway name
         x = "",  # Remove x-axis label
         y = "GSVA Score") +
    theme(text = element_text(size = 14, family = "Arial"),  # Ensure all text uses Arial
          plot.title = element_text(hjust = 0.5, size = 14),  # Center title
          axis.text.x = element_text(size = 16, color = "black"),  # Increase x-axis text size and set color to black
          axis.text.y = element_text(size = 14, color = "black"),  # Set y-axis text color to black
          axis.title.x = element_text(size = 14, color = "black"),  # Set x-axis title color to black
          axis.title.y = element_text(size = 14, color = "black"),  # Set y-axis title color to black
          legend.position = "none")  
  
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), p,
         width = fig.width, height = fig.height, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), p,
         width = fig.width, height = fig.height, units = "in")
  
  return(p)
  
}



plot_gsva_boxplot_mutil_5<-  function(gsva_matrix, condition_list_label, pathway_name,
                                      figure_folder, file_name, 
                                      fig.height = 6, fig.width = 8,
                                      reference_group_1, compare_group_1,
                                      reference_group_2, compare_group_2,
                                      reference_group_3, compare_group_3,
                                      reference_group_4, compare_group_4,
                                      reference_group_5, compare_group_5){
  
  # Sample list sub
  sample_info <- condition_list_label %>%
    mutate(sample = rownames(.)) 
  
  # Convert gsva_matrix to a data frame and reshape
  plot_df <- as.data.frame(gsva_matrix) %>%
    rownames_to_column(var = "pathway") %>%
    filter(pathway == pathway_name) %>%  # Select the specific pathway
    pivot_longer(cols = -pathway, names_to = "sample", values_to = "GSVA_score") %>%
    dplyr::select(-pathway)  %>% 
    filter(sample %in% sample_info$sample)  %>%
    left_join(sample_info, by = "sample") 
  
  # Ensure colors are mapped to exact group names
  # color_palette <- setNames(c("#10d5da", "#fe867f"), c(reference_group, compare_group))
  
  # Process pathway name: remove first part, capitalize first letter, replace underscores with spaces
  formatted_title <- pathway_name %>%
    str_remove("^[^_]+_") %>%  # Remove everything before the first underscore
    str_to_lower() %>%  # Convert everything to lowercase
    str_replace_all("_", " ") %>%  # Replace underscores with spaces
    str_to_sentence()  # Capitalize first letter
  
  # sub the rna with RNA in formatted_title
  formatted_title <- str_replace_all(formatted_title, "rna", "RNA")
  
  # Create the box plot with scatter overlay
  p<-ggplot(plot_df, aes(x = group, y = GSVA_score)) +
    geom_boxplot(aes(fill = group), alpha = 0.9, outlier.shape = NA, color = "black") +  # Box plot with fill color
    geom_jitter(aes(fill = group), shape = 21, width = 0.2, size = 3, alpha = 0.9, color = "black") +  # Scatter points with fill color
    # stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +  # Mean point
    geom_signif(comparisons = list(c(reference_group_1, compare_group_1),
                                   c(reference_group_2, compare_group_2),
                                   c(reference_group_3, compare_group_3),
                                   c(reference_group_4, compare_group_4),
                                   c(reference_group_5, compare_group_5)),
                test = "t.test", 
                map_signif_level = TRUE) +  # Add significance annotation
    # scale_fill_manual(values = color_palette) +  # Apply custom colors to boxes & scatter dots
    theme_minimal(base_family = "Arial") +  # Use Arial font for all text
    labs(title = formatted_title,  # Use formatted pathway name
         x = "",  # Remove x-axis label
         y = "GSVA Score") +
    theme(text = element_text(size = 14, family = "Arial"),  # Ensure all text uses Arial
          plot.title = element_text(hjust = 0.5, size = 14),  # Center title
          axis.text.x = element_text(size = 16, color = "black"),  # Increase x-axis text size and set color to black
          axis.text.y = element_text(size = 14, color = "black"),  # Set y-axis text color to black
          axis.title.x = element_text(size = 14, color = "black"),  # Set x-axis title color to black
          axis.title.y = element_text(size = 14, color = "black"),  # Set y-axis title color to black
          legend.position = "none")  
  
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), p,
         width = fig.width, height = fig.height, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), p,
         width = fig.width, height = fig.height, units = "in")
  
  return(p)
  
}



