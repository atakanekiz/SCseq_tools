CIPR <- function(input_dat, 
                 comp_method = "dot_product",
                 reference = "immgen",
                 custom_ref_dat_path = NULL,
                 custom_ref_annot_path = NULL,
                 plot_ind = F, 
                 plot_top = T, 
                 top_num = 5, 
                 save_png = F, 
                 global_plot_obj = T,
                 global_results_obj = T){
  
  # adapted from v4
  
  suppressMessages({
    require(ggpubr)
    require(gtools)
    require(tibble)
    require(dplyr)
  })
  
  
  
  # Define column names to allow flexibility in case and close matches in column names
  gene_column <- grep("gene", colnames(input_dat), ignore.case = T, value = T)
  logFC_column <- grep("logfc", colnames(input_dat), ignore.case = T, value = T)
  cluster_column <- grep("cluster", colnames(input_dat), ignore.case = T, value = T)
  
  
  # Convert gene symbols to lower case letters to allow mouse-vs-human comparisons
  input_dat[,gene_column] <- tolower(input_dat[,gene_column])
  
  if(reference == "immgen"){
    
    message("Reading ImmGen (v1+v2) reference data")
    
    # Read reference dataset
    ref_log <- readRDS(url("https://github.com/atakanekiz/CIPR_development/blob/master/CIPR_v4/data/immgen_recalc_ratio.rds?raw=true"))
    
    
    # Read immgen annotation file for explanations of cell types
    ref_annot <- readRDS(url("https://github.com/atakanekiz/CIPR_development/blob/master/CIPR_v4/data/immgen_annot.rds?raw=true"))
    
    
  } else {
    
    message("Reading custom reference data")
    
    # Read reference dataset
    ref_log <- readRDS(custom_ref_dat_path)
    
    
    # Read immgen annotation file for explanations of cell types
    ref_annot <- readRDS(custom_ref_annot_path)
    
  }
  
  ref_gene_column <- grep("gene", colnames(ref_log), ignore.case = T, value = T)
  
  ################################################################################################################################
  # Define a reactive cluster object that will store cluster information
  # Get the clusters and sort them in incrementing order from cluster column
  # This is needed to generate results per cluster
  
  clusters <- gtools::mixedsort(
    levels(
      as.factor(
        pull(input_dat, grep("cluster", x = colnames(input_dat), ignore.case = T, value = T)
        )
      )
    )
  )
  
  
  
  
  
  
  
  ################################################################################################################################
  # Compare user_data against reference file
  
  message("Analyzing cluster signatures")
  
  
  if(comp_method == "dot_product"){
    
    # Initiate a master data frame to store the results
    master_df <- data.frame()
    
    # Iterate over clusters to calculate a distinct identity score for each reference cell type
    for (i in clusters) {
      
      # Increment the progress bar, and update the detail text.
      message(paste("Analyzing cluster", i))
      
      # Subset on the cluster in iteration
      sel_clst <- input_dat %>%
        filter(!!rlang::sym(cluster_column) == i) %>%
        select(c(!!sym(gene_column), !!sym(logFC_column)))
      
      
      # Merge SCseq cluster log FC value with immgen log FC for shared genes
      merged <- merge(sel_clst, ref_log, by.x = gene_column, by.y = ref_gene_column)
      
      
      # Calculate a scoring matrix by multiplying log changes of clusters and immgen cells
      reference_scoring <- data.frame(apply(merged[,3:dim(merged)[2]],2,function(x){x*merged[,2]}), check.names = FALSE)
      
      # Calculate the aggregate score of each immgen cell type by adding
      score_sum <- colSums(reference_scoring)
      
      # Store identity scores in a data frame
      df <- data.frame(identity_score = score_sum)
      
      df <- rownames_to_column(df, var="reference_id")
      
      
      df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))
      
      
      # Store cluster information in a column
      df$cluster <- i
      
      # Add confidence-of-prediction calculations here and append to the df
      # Calculate the mean and standard deviation of the aggregate scores per reference cell type
      mean_score_sum <- mean(df$identity_score)
      score_sum_sd <- sd(df$identity_score)
      
      # Calculate the distance of the identity score from population mean (how many std devs apart?)
      df$z_score <- (df$identity_score - mean_score_sum)/score_sum_sd
      
      # Calculate the proportion of the genes changing in the same direction between unknown cluster and reference cell type
      df$percent_pos_correlation <- {
        
        ngenes <- dim(reference_scoring)[1]
        
        pos_corr_vector <- numeric()
        
        for(i in 1:dim(reference_scoring)[2]){
          
          # Calculate number of genes positively correlated (upregulated or downregulated in both unk cluster and reference)
          pos_cor <- ( sum(reference_scoring[, i] > 0) / ngenes ) * 100
          
          pos_corr_vector <- c(pos_corr_vector, pos_cor)
          
        } #close for loop
        
        pos_corr_vector
        
      } # close expression 
      
      
      # Add calculation results under the master data frame to have a composite results file
      master_df <- rbind(master_df,df)
      
      
      
    } # close for loop that iterates over clusters
    
  } else {  ################### Correlation methods ###########################################################

    # Initiate master data frame to store results
    master_df <- data.frame()


    # Iterate analysis for each cluster. The loop below will calculate a distinct correlation
    # coefficient for each cluster-reference cell pairs
    for (i in clusters) {


      trim_dat <- input_dat %>%
        filter(!!rlang::sym(cluster_column) == i)

      dat_genes <- trim_dat[gene_column] %>% pull() %>% as.character
      ref_genes <- ref_log[ref_gene_column] %>% pull() %>% as.character

      common_genes <- intersect(dat_genes, ref_genes)


      trim_dat <- trim_dat %>%
        filter(!!rlang::sym(gene_column) %in% common_genes) %>%
        arrange(!!rlang::sym(gene_column)) %>%
        select(- !!rlang::sym(gene_column))


      trim_ref <- ref_log %>%
        filter(!!rlang::sym(ref_gene_column) %in% common_genes) %>%
        arrange(!!rlang::sym(ref_gene_column)) %>%
        select(- !!rlang::sym(ref_gene_column))


      # Calculate correlation between the the cluster (single column in trimmed input data) and each of the
      # reference cell subsets (columns of the trimmed reference data)
      cor_df <- cor(trim_dat[logFC_column], trim_ref, method = comp_method)

      # Store results in a data frame
      df <- data.frame(identity_score = cor_df[1,])

      df <- rownames_to_column(df, var="reference_id")

      # Combine results with reference annotations
      if(input$sel_reference == "ImmGen"){

        df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))





      } else if (input$sel_reference == "Custom" & !is.null(input$annot_file)){

        df <- left_join(df, ref_annot, by=c("reference_id" = "short_name"))


      } else if(input$sel_reference == "Custom" & is.null(input$annot_file)){

        # Fill in with reminder if annotation file is not updated
        df$reference_cell_type <- rep("Upload annotation file", dim(ref_log)[2]-1)
        df$short_name <- colnames(ref_log)[!colnames(ref_log) %in% ref_gene_column]
        df$long_name <- rep("Upload annotation file", dim(ref_log)[2]-1)
        df$description <- rep("Upload annotation file", dim(ref_log)[2]-1)

      }



      # Store cluster information in a column
      df$cluster <- i

      # Add confidence-of-prediction calculations here and append to the df
      # Calculate the mean and standard deviation of the aggregate scores per reference cell type
      mean_cor_coeff <- mean(df$identity_score)
      cor_coeff_sd <- sd(df$identity_score)

      # Calculate the distance of the identity score from population mean (how many std devs apart?)
      df$z_score <- (df$identity_score - mean_cor_coeff)/cor_coeff_sd

      # Add all the results to the master data frame
      master_df <- rbind(master_df, df)


    } # close for loop that iterates over clusters

}
    if(global_results_obj == T) CIPR_results <<- master_df
    
  #prep individual plots
  if(plot_ind == T){
    
    ind_clu_plots <- list()      
    
    for (i in clusters) {
      
      
      # Extract results calculated for individual clusters
      df_plot <- master_df %>%
        filter(cluster == i)
      
      # Calculate mean and sd deviation for adding confidence bands to graphs
      score_mean <- mean(df_plot$identity_score)
      score_sd <- sd(df_plot$identity_score)
      
      
      plotname <- paste("cluster", i, sep="")
      
      # Plot identity scores per cluster per reference cell type and add confidence bands
      ind_clu_plots[[plotname]] <- ggdotplot(df_plot, x = "reference_id", y="identity_score", 
                                             fill = "reference_cell_type", xlab=F, ylab="Reference identity score",
                                             font.y = c(14, "bold", "black"), size=1, x.text.angle=90,
                                             title = paste("Cluster:",my_i), font.title = c(15, "bold.italic"),
                                             font.legend = c(15, "plain", "black"))+
        theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=1))+
        geom_hline(yintercept=score_mean)+
        annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                 ymin = score_mean-score_sd, ymax = score_mean+score_sd,
                 fill = "gray50", alpha = .1)+
        annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                 ymin = score_mean-2*score_sd, ymax = score_mean+2*score_sd,
                 fill = "gray50", alpha = .1)
      
      
    }
    
    if(global_plot_obj == T) ind_clu_plots <<-ind_clu_plots
    
    
    if(save_png == T) {
      ggexport(filename = "CIPR_individual_clusters.png", plotlist = ind_clu_plots, ncol = 1, width = 1800, height = 360 * length(clusters))
    }
    else {
      print(ggarrange(plotlist = ind_clu_plots, ncol = 1, common.legend = T))
    }
    
  }
  
  ################################################################################################################################
  # Prepare top5 summary plots
  # This plot will show the 5 highest scoring reference cell types for each cluster.
  
  if(plot_top == T){
    
    # Extract top5 hits from the reuslts
    top_df <- master_df %>%
      group_by(cluster) %>%    #cluster
      top_n(top_num, wt = identity_score) %>%
      arrange(as.numeric(cluster), desc(identity_score))
    
    # Index variable helps keeping the results for clusters separate and helps ordered outputs
    top_df$index <- 1:nrow(top_df)
    
    
    # Order clusters levels for ordered plotting
    ordered_cluster_levels <- gtools::mixedsort(levels(as.factor(top_df$cluster)))
    
    
    top_df$cluster <- factor(top_df$cluster, levels = ordered_cluster_levels)
    
    
    
    # Extract relevant columns 
    top_df <- select(top_df, cluster,
                     reference_cell_type,
                     reference_id,
                     long_name,
                     description,
                     identity_score,
                     index, everything())
    
    if(global_results_obj == T) CIPR_top_df <<- top_df
    
    p <- ggdotplot(top_df, x="index", y="identity_score", 
                   fill = "cluster", size=1, x.text.angle=90, 
                   font.legend = c(15, "plain", "black")) +
      scale_x_discrete(labels=top_df$reference_id)+
      theme(axis.text.x = element_text(vjust=0.5, hjust=1))
    
    if(global_plot_obj == T) top_plots <<- p
    
    if(save_png == T) {
      
      ggexport(p, filename = "CIPR_top_hits.png", ncol = 1, width = 150 * length(clusters), height = 300)
      
    } else {
      
      print(p)
      
    }
    
  }
  
} # close function
