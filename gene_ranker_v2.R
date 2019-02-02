# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# A function to rank genes in a gene expression dataframe
# Uses dplyr
###################################################################################################################################################################
gene_ranker <- function(exprs = NULL, # Expression data frame (rows are cells, columns are genes and metadata)
                        pos_marker = NULL, # Genes to positively gate cells (cells expressing these markers will be subsetted)
                        neg_marker = NULL, # Genes to negatively gate cells (cells expressing these markers will be discarded)
                        sample_id, # Which cells will be selected as 'sample' (ie, the direction of rank ordering. Regex based string recognition
                        sample_cluster = NULL, # Cell clusters to include in analysis for the sample subset
                        reference_id, # Which cells will be selected as 'reference' (ie, the direction of rank ordering)
                        reference_cluster = NULL, # Cell clusters to include in analysis for the reference subset
                        method = "s2n", # can also take the following values ("ttest", "difference", "ratio", "welch", "mwt", "bws"), # Method for ranking the genes. For MWT, see PMID: 18344518  
                        save_ranked_object = F # Set it to TRUE if you want to assign a global variable for ranking data
                        ){
  
  require(dplyr)
  require(stringr)
  require(tibble)
  
  
  # Subset cells that express markers of interest at any level (>0)
  if(!is.null(pos_marker)) {exprs <- filter_(exprs, paste(pos_marker , "!= 0", collapse = " & "))}
  
  # Discard cells that express a marker gene which we want to negatively gate in our analyses
  if(!is.null(neg_marker)) {exprs <- filter_(exprs, paste(neg_marker , "== 0", collapse = " & "))}
  
  # Subset cells to serve as sample in comparison (sample vs ref comparison i.e. left side of the GSEA plots later on)
  sample_df <- filter(exprs, Sample %in% sample_id)
  

  
  # Subset clusters of interest
  if(!is.null(sample_cluster)) {
    
    sample_df <- filter(sample_df, Cluster %in% sample_cluster)
    
    }
  
  
  
  
  # Subset cells to serve as reference in comparison (sample vs ref comparison i.e. right side of the GSEA plots later on)
  reference_df <- filter(exprs, Sample %in% reference_id)
  
  

  
  # Subset clusters of interest
  if(!is.null(reference_cluster)) {

    reference_df <- filter(reference_df, Cluster %in% reference_cluster)
    
    }
  

  if(method == "mwt"){
    
    require(mwt)
    require(BiocGenerics)
    
    discard_from_sample <- as.logical(colSums(sample_df == 0, na.rm = T))
    discard_from_ref <- as.logical(colSums(reference_df == 0, na.rm = T))
    discard_from_both <- discard_from_sample | discard_from_ref
    sample_df <- sample_df[,!discard_from_both]
    reference_df <- reference_df[,!discard_from_both]
   
    sample_df <- add_column(sample_df, Comparison = "Sample", .after = 0)
    
    reference_df <-add_column(reference_df, Comparison="Reference", .after=0)
    
    
    
    
    trimmed_df <- rbind(sample_df, reference_df)
    group_colnames <- as.factor(trimmed_df$Comparison)
    trimmed_df <- trimmed_df %>%
      dplyr::select(-Sample, -Cluster, -Cell_id, -Comparison)%>%
      t()
    
    colnames(trimmed_df) <- group_colnames
    
    
  
    mwt_results <- mwt(trimmed_df, grp = group_colnames)
    
  } else {
    
    
    # Report which cells are being analyzed
    print(paste(dim(sample_df)[1], "cells with the following annotations are included in the 'sample' group (left side of GSEA plot)"))
    print(paste("Sample:",paste(levels(droplevels(as.factor(sample_df$Sample))), collapse = ", ")))
    print(paste("Sample clusters in analysis:", paste(levels(droplevels(as.factor(sample_df$Cluster))), collapse = ", ")))
    
    
    print(paste(dim(reference_df)[1], "cells with the following annotations are included in the 'reference' group (right side of GSEA plot)"))
    print(paste("Reference:", paste(levels(droplevels(as.factor(reference_df$Sample))), collapse = ", ")))
    print(paste("Reference clusters in analysis:", paste(levels(droplevels(as.factor(reference_df$Cluster))), collapse = ", ")))
    
    
  
    # Get rid of unnecessary columns
    sample_df <- sample_df[,!colnames(sample_df) %in% c("Sample", "Cluster", "Cell_id")]
    reference_df <- reference_df[,!colnames(reference_df) %in% c("Sample", "Cluster", "Cell_id")]
    
    sample_mean <- colMeans(sample_df)
    sample_sd <- apply(sample_df, MARGIN = 2, FUN = sd)
    sample_var <- apply(sample_df, MARGIN = 2, FUN = var)
    sample_n <- dim(sample_df)[1]
    
    reference_mean <- colMeans(reference_df)
    reference_sd <- apply(reference_df, MARGIN = 2, FUN = sd)
    reference_var <- apply(reference_df, MARGIN = 2, FUN = var)
    reference_n <- dim(reference_df)[1]
    
  }
  
  

  
  
  

  
  

  
  
  
  if(method == "s2n"){
    
    ranking <- (sample_mean - reference_mean) / (sample_sd + reference_sd)
    
  } else if(method == "ttest"){
    
    ranking <- (sample_mean - reference_mean) / sqrt( ((sample_sd^2)/sample_n) + ((reference_sd^2)/reference_n) )
    
  } else if(method == "difference") {
    
    ranking <- sample_mean - reference_mean
    
  } else if(method == "ratio") {
    
    ranking <- sample_mean / reference_mean
    
  } else if (method == "welch") {
    
    ranking <- (sample_mean - reference_mean) / sqrt( ((sample_var^2)/sample_n) + ((reference_var^2)/reference_n) ) 
    
  } else if(method == "mwt"){
    
    ranking <- mwt_results$MWT
    # names(ranking) <- rownames(trimmed_df) #not needed?
    
  } else if (method == "bws"){
    
    require(BWStest)
    require(dplyr)
    
    if(colnames(sample_df) != colnames(reference_df)){
      
      print("Different column (gene) names between sample and reference. Using shared genes for ranking")
      
      shared <- intersect(colnames(sample_df), colnames(reference_df))
      
      sample_df <- sample_df[,shared]
      reference_df <- reference_df[,shared]
      
      #ensure match
      if(colnames(sample_df) != colnames(reference_df)) error("Column (gene) names still do not match between sample and reference")
      
      }
    
    ranking <- numeric() 
    
    for(i in colnames(sample_df)){
      
      sample_gene <- pull(sample_df, i)
      reference_gene <- pull(reference_df,i)
      
      ranking[i] <- bws_stat(sample_gene, reference_gene)
      
    }
    
  } else {stop('Select one of the following statistical tests: s2n, ttest, difference, ratio, welch, mwt, bws')}
  
  
  ranking <- subset(ranking, c(!is.na(ranking) & !is.infinite(ranking)))
  ranking <- sort(ranking, decreasing = T)
  
  if(save_ranked_object == T) {assign("ranked_geneset", ranking, envir = .GlobalEnv)} # Useful for iterative analysis when executed in pipeline
  
  ranking # Return ranked gene set
    
}





# # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test # Test
# 

# exprs <- readRDS("aging_exprs.rds")
# 
# gene_ranker(exprs, 
#              sample_id = "Young \\(WT\\)", 
#              reference_id = "Aged \\(WT\\)",
#              sample_cluster = "Act Cd8+",
#              reference_cluster = "Act Cd8+")





# s2n_ranked <- gene_ranker(exprs = exprs,
#                           pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                           sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                           method = "s2n")
# 
# ttest_ranked <- gene_ranker(exprs = exprs, 
#                             pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                             sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                             method = "ttest")
# 
# welch_ranked <- gene_ranker(exprs = exprs, 
#                             pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                             sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                             method = "welch")
# 
# difference_ranked <- gene_ranker(exprs = exprs, 
#                                  pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                                  sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                                  method = "difference")
# 
# ratio_ranked <- gene_ranker(exprs = exprs, 
#                             pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                             sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                             method = "ratio")
# 
# mwt_ranked <- gene_ranker(exprs = exprs, 
#                           pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                           sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                           method = "mwt")
# 
# bws_ranked <- gene_ranker(exprs = exprs, 
#                           pos_marker = c("Cd3e", "Cd8a"), neg_marker = "Cd4",
#                           sample_id = "Young \\(WT\\)", reference_id = "Aged \\(WT\\)",
#                           method = "bws")
# 
# 
# # Compare results
# 
# res_list<-list(s2n_ranked = s2n_ranked, 
#                ttest_ranked = ttest_ranked,
#                welch_ranked = welch_ranked, 
#                difference_ranked = difference_ranked, 
#                ratio_ranked = ratio_ranked,
#                bws_ranked=bws_ranked) 
#                # mwt_ranked = mwt_ranked)
# 
# common_genes <- Reduce(intersect, list(names(s2n_ranked),
#                                        names(ttest_ranked),
#                                        names(welch_ranked),
#                                        names(difference_ranked),
#                                        names(ratio_ranked),
#                                        names(bws_ranked)))
# 
# trimmed <- lapply(res_list, FUN=function(x) x[common_genes])
# 
# 
# combined_res <- do.call(cbind.data.frame, trimmed)
# 
# 
# library(corrplot)
# library(RColorBrewer)
# M <-cor(combined_res)
# corrplot(M, type="upper", order="hclust",
#          col=brewer.pal(n=10, name="RdYlBu"), method = "pie")
