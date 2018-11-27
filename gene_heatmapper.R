# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# exprs <- readRDS("aging_exprs.rds")


gene_heatmapper <- function(exprs, # expression dataframe (generated with df_extractor function)
                            genes_to_plot, # genes to plot 
                            clusters_to_plot = NULL, # select clusters to be subsetted for plot
                            pos_marker = NULL, # select cells expressing these markers. Provide a comma separated string of gene names 
                            neg_marker = NULL, # select cells not expressing these markers. Provide a comma separated string of gene names 
                            annotate = c("Sample", "Cluster"),
                            scale = "none",
                            expand_scale = F,
                            cluster_rows = T,
                            cluster_cols = T,
                            clustering_distance_rows = "euclidean",
                            clustering_distance_cols = "euclidean",
                            clustering_method = "complete",
                            save_pdf = T, # Save results in a pdf file in the workspace.
                            append_to_filename = "",
                            output_plot = F, # Set to TRUE if you'd like to show graph in RMD or viewer
                            ... # Not used currently, but can be included to pass other arguments to graphing functions
)


{
  
  require(pheatmap)
  require(dplyr)
  require(stringr)
  require(RColorBrewer)
  # require(purrr)
  # require(colorspace)
  require(gtools)
  
  colnames(exprs) <- str_replace_all(colnames(exprs), "[:punct:]|[:space:]", "_")
  genes_to_plot <- str_replace_all(genes_to_plot, "[:punct:]|[:space:]", "_")
  genes_to_plot <- genes_to_plot[genes_to_plot %in% colnames(exprs)]
  genes_to_plot <- mixedsort(genes_to_plot)
  
  # Subset clusters of interest
  if(!is.null(clusters_to_plot)) {exprs <- filter(exprs, str_detect(Cluster, regex(clusters_to_plot, ignore_case=T, comments = F)))}
  
  # Subset cells that express markers of interest at any level (>0)
  if(!is.null(pos_marker)) {
    
    pos_marker_split <- trimws(strsplit(pos_marker, split = ",")[[1]])
    exprs <- filter_(exprs, paste(pos_marker_split , "!= 0", collapse = " & "))}
  
  # Discard cells that express a marker gene which we want to negatively gate in our analyses
  if(!is.null(neg_marker)) {
    neg_marker_split <- trimws(strsplit(neg_marker, split = ",")[[1]])
    exprs <- filter_(exprs, paste(neg_marker_split , "== 0", collapse = " & "))}
  
  
  # Report which cells are being analyzed
  message(paste(dim(exprs)[1], "cells with the following annotations will be plotted"))
  message("\nSample:\n", 
          paste(levels(droplevels(as.factor(exprs$Sample))), collapse = ", "), "\n",
          "\nSample clusters in analysis:\n", 
          paste(levels(droplevels(as.factor(exprs$Cluster))), collapse = ", "), "\n\n")
  
  
  omitted_genes <- character()
  
  for (gene in genes_to_plot){
    
    if(sum(exprs[,gene]) == 0){
      
      omitted_genes <- c(omitted_genes, gene)
      
    }
    
  }
  
  genes_to_plot <- genes_to_plot[!genes_to_plot %in% omitted_genes]
  
  if(length(genes_to_plot) == 0) stop("Requested genes have no expression in the subsetted data frame")
  
  if(length(omitted_genes) != 0) message(paste("Following genes are omitted due to zero expression: ", paste(omitted_genes, collapse = ", ")))
  

  exprs_mtx <- as.matrix(select(exprs, genes_to_plot))
  
  if(!is.na(annotate)){  
    exprs_metadata <- as.data.frame(select(exprs, annotate)) 
    } else {exprs_metadata = NA}
    
  # Plot heatmap
  
  if(expand_scale == T) {
    
    stat_summary <- summary(as.numeric(as.matrix(exprs)))
    
    rounded_mean <- ceiling(as.numeric(stat_summary["Mean"]))
    min_value <- as.numeric(stat_summary["Min."])
    max_value <- as.numeric(stat_summary["Max."])
    
    breaks = seq(min_value, rounded_mean ,by=(rounded_mean-min_value)/100)
    
    #Add outliers
    breaks = append(breaks, max_value)
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(length(breaks))
      
  } else {
    
    breaks = NA
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
    
    }
  
  
  if(output_plot == T) {
    
    pheatmap(exprs_mtx, annotation_row = exprs_metadata,
           breaks = breaks, color=color,
           scale = scale, 
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           clustering_distance_rows = clustering_distance_rows,
           clustering_distance_cols = clustering_distance_cols,
           clustering_method = clustering_method,
           border_color = NA)
    
  }
  
  if(save_pdf == T){
    
    arg_list <- list(cluster = clusters_to_plot, pos = pos_marker, neg = neg_marker)
    select_non_null <- !sapply(arg_list, is.null)
    
    if(sum(select_non_null) == 0) {
      
      filename <- paste0("Unsubsetted_data_plots__", stat_method, ".pdf")
      
    } else {
      
      filename <- paste(arg_list[select_non_null], names(arg_list[select_non_null]), collapse = "__")
      filename <- paste0(filename,"__", clustering_method, append_to_filename, ".pdf")
    }
    
    
    
    pheatmap(exprs_mtx, annotation_row = exprs_metadata,
             breaks = breaks, color=color,
             scale = scale, 
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             clustering_distance_rows = clustering_distance_rows,
             clustering_distance_cols = clustering_distance_cols,
             clustering_method = clustering_method,
             filename = filename,
             border_color = NA)
    
  } 
  
}



############################################### TEST ################################################# 

ms_immune <- read.delim("mouse_immune_process.txt", header = F)
ms_immune <- ms_immune$V1

gene_heatmapper(exprs, annotate = NA,
                genes_to_plot = ms_immune,
                clusters_to_plot = "Act Cd8+", 
                save_pdf = T,
                expand_scale = T,
                cluster_cols = T,
                cluster_rows = T)

