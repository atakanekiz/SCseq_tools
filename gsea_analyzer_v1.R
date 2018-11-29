library(fgsea)
library(dplyr)
library(grid)
library(ggplot2)

# Read molecular signatures database (MSigDB) gene lists. ALL list contains the other ones plus extra lists
msig_ALL <- gmtPathways("msigdb_genelists/msigdb.v6.2.symbols.gmt")
msig_CURATED <- gmtPathways("msigdb_genelists/c2.all.v6.2.symbols.gmt")
msig_MOTIF <- gmtPathways("msigdb_genelists/c3.all.v6.2.symbols.gmt")
msig_GO <- gmtPathways("msigdb_genelists/c5.all.v6.2.symbols.gmt")
msig_IMMUNE <- gmtPathways("msigdb_genelists/c7.all.v6.2.symbols.gmt")
msig_HALLMARK <- gmtPathways("msigdb_genelists/h.all.v6.2.symbols.gmt")

########################################################################################################################
########################################################################################################################
# Create function to extract rank gene expression between two groups
# This method uses S2N metric (mean(a)-mean(b)/sd(a)+sd(b))
# S2N metric as written above result ranks 'a' in comparison to 'b' (left side of the plots)
# To use the function, you can write marker information and/or sample vs reference information
# Genes are ranked higher in the sample input. Therefore, 'sample' will appear on the left side 
#  of GSEA columns (reference is on the right side)
rank_generator <- function(exprs = NULL, marker1 = NULL, marker2 = NULL,
                           samp_id1 = NULL,samp_id2 = NULL,
                           ref_id1 = NULL,ref_id2 = NULL) {
  
  if (!is.null(marker1)) {
    
    marker1_pos <- which(exprs[marker1, ] != 0)
    marker_col_ID <- marker1_pos
    
    if (!is.null(marker2)) {
      
      marker2_pos <- which(exprs[marker2, ] != 0)
      marker_col_ID <- intersect(marker_col_ID, marker2_pos)
    }
    
    
    samp_1 <- which(grepl(samp_id1,  colnames(exprs), ignore.case = T))
    if (!is.null(samp_id2)) {
      samp_2 <- which(grepl(samp_id2,  colnames(exprs), ignore.case = T))
      samp <- intersect(samp_1, samp_2, marker_col_ID)
      message(paste(length(samp), "cells with the following annotations are included in the 'sample' group (left side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,samp]))))))
    } else {
      samp <-intersect(samp_1, marker_col_ID)
      message(paste(length(samp), "cells with the following annotations are included in the 'sample' group (left side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,samp]))))))
    }
    
    ref_1 <- which(grepl(ref_id1,  colnames(exprs), ignore.case = T))
    if(!is.null(ref_id2)){
      ref_2 <- which(grepl(ref_id2,  colnames(exprs), ignore.case = T))
      ref <- intersect(ref_1, ref_2, marker_col_ID)
      message(paste(length(ref), "cells with the following annotations are included in the 'reference' group (right side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,ref]))))))
    }  else {
      ref<- intersect(ref_1, marker_col_ID)
      message(paste(length(ref), "cells with the following annotations are included in the 'reference' group (right side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,ref]))))))
    }
    
  } else {
    
    samp_1 <- which(grepl(samp_id1,  colnames(exprs), ignore.case = T))
    if (!is.null(samp_id2)) {
      samp_2 <- which(grepl(samp_id2,  colnames(exprs), ignore.case = T))
      samp <- intersect(samp_1, samp_2)
      message(paste(length(samp), "cells with the following annotations are included in the 'sample' group (left side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,samp]))))))
    } else {
      samp <- samp_1
      message(paste(length(samp), "cells with the following annotations are included in the 'sample' group (left side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,samp]))))))
    }
    
    ref_1 <- which(grepl(ref_id1,  colnames(exprs), ignore.case = T))
    if(!is.null(ref_id2)){
      ref_2 <- which(grepl(ref_id2,  colnames(exprs), ignore.case = T))
      ref <- intersect(ref_1, ref_2)
      message(paste(length(ref), "cells with the following annotations are included in the 'reference' group (right side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,ref]))))))
    } else { 
      ref <- ref_1
      message(paste(length(ref), "cells with the following annotations are included in the 'reference' group (right side of GSEA plot)", sep = " "))
      print(levels(as.factor(sub("_[^_]+$", "", (colnames(exprs[,ref]))))))
    }
  }
  
  assign("sample_group", value = samp, envir = .GlobalEnv)
  assign("reference_group", value = ref, envir = .GlobalEnv)
  
  ref_avg <- rowMeans(exprs[,ref])
  ref_stdev <- apply(exprs[,ref], 1, sd)
  samp_avg <- rowMeans(exprs[,samp])
  samp_stdev <- apply(exprs[,samp], 1, sd)
  
  ranking <- (samp_avg - ref_avg)/(samp_stdev + ref_stdev)
  ranking <- subset(ranking, c(!is.na(ranking) & !is.infinite(ranking)))
  ranking <- sort(ranking, decreasing = T)
  ranking
}



##############################################################################################################
# # Example usage of rank_generator function
# f480_d12_rnk <- rank_generator(exprs = exprs, marker1 = "ADGRE1", samp_id1 ="d12wt", ref_id1 = "d12ko")
# 
# act_cd8_d12_rnk <- rank_generator(exprs = exprs, 
#                                   samp_id1 = "activated_cd8", samp_id2 = "d12wt",
#                                   ref_id1 = "activated_cd8", ref_id2 = "d12ko")
# nai_cd8_d12_rnk <- rank_generator(exprs = exprs, 
#                                   samp_id1 = "naive_cd8", samp_id2 = "d12wt",
#                                   ref_id1 = "naive_cd8", ref_id2 = "d12ko")
# 
# cd3cd8_d12_WTvsKO_rnk <- rank_generator(exprs = exprs, marker1 = "CD8A", marker2 = "CD3E", samp_id1 ="d12wt", ref_id1 = "d12ko")
# 
# 
# # Feed the ranked genes into enrichment function 
# f480_d12_gsea <- fgsea(pathways = msig_HALLMARK, stats = f480_d12_rnk, 
#                        nperm = 10000, minSize = 50, maxSize = 500)
# act_cd8_d12_gsea <- fgsea(pathways = msig_HALLMARK, stats = act_cd8_d12_rnk, 
#                           nperm = 10000, minSize = 50, maxSize = 500)
# nai_cd8_d12_gsea <- fgsea(pathways = msig_HALLMARK, stats = nai_cd8_d12_rnk, 
#                           nperm = 10000, minSize = 50, maxSize = 500)
# cd3cd8_d12_WTvsKO_gsea <- fgsea(pathways = msig_HALLMARK, stats = cd3cd8_d12_WTvsKO_rnk, 
#                                 nperm = 10000, minSize = 50, maxSize = 500) 
# 
# 
# # Rank lists were generated based on sample/ref ratio. Hence, the left side of the enrichment
# # plots (low ranked genes) represent WT group.
# plotEnrichment(msig_HALLMARK[[grep("INTERFERON_GAMMA_RESPONSE", names(msig_HALLMARK), value = T)]],
#                f480_d12_rnk) + labs(title=grep("INTERFERON_GAMMA_RESPONSE", names(msig_HALLMARK), value = T))
##############################################################################################################
##############################################################################################################








##############################################################################################################
##############################################################################################################
# Create top_plotter function to speed up plotting  
# Use top plotter function with previously generated "gene rankings" and "gsea results"
top_plotter <- function(gsea_results = NULL, ranked_genes = NULL, gene_set = NULL, 
                        top_n = 10, gseaParam = 1) {
  topPathwaysUp <- gsea_results[ES > 0][head(order(pval), n=top_n), pathway]
  topPathwaysDown <- gsea_results[ES < 0][head(order(pval), n=top_n), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(gene_set[topPathways], ranked_genes, gsea_results, 
                gseaParam = gseaParam)
  
}

##############################################################################################################
# # Example usage of top_plotter function
# top_plotter(gsea_results = f480_d12_gsea, ranked_genes = f480_d12_rnk, gene_set = msig_HALLMARK, top_n = 5, gseaParam = 0.4)
# top_plotter(gsea_results = act_cd8_d12_gsea, ranked_genes = act_cd8_d12_rnk, gene_set = msig_HALLMARK, top_n = 5, gseaParam = 0.4)
##############################################################################################################
##############################################################################################################



##############################################################################################################
##############################################################################################################
################################# The ultimate automated function ############################################

# Define master_gsea function that does the following tasks automatically:
# Extract ranked genes in sample vs reference (S2N metric) (remember, sample will appear on the left of GSE plot)
# Perform fgsea analysis with desired gene list (gene_set parameter)
# Report top hits of the enrichment results (the number of plotted graphs can be adjusted)
# Alternatively, you can plot an individual enrichment plot with annotations showing p-value and NES values
# Length of the gene hits (bars) can be adjusted with gseaParam function (pick 0.2-0.5 for >10-15 graphs)
# This function will also ask the user which pathway needs to be plotted in the case of multiple string matches.
# Parameters are passed as character class. Both lower and upper case values work.
# The function also assigns global environment variables (ranked_gene, gene lists used, 
# gsea_res, and some other intermediate variables to sort out the multiple match situation)

master_gsea <- function(exprs = NULL, marker1 = NULL, marker2 = NULL,
                        samp_id1 = NULL, samp_id2 = NULL,
                        ref_id1 = NULL, ref_id2 = NULL, 
                        gene_set = NULL, nperm = 10000, minSize = 50, maxSize = 500,
                        top_n = 10, gseaParam = 1, plot_individual = NULL, append_title = F){
  
  ranked_genes <- rank_generator(exprs = exprs, marker1 = marker1, marker2 = marker2,
                                 samp_id1 = samp_id1, samp_id2 = samp_id2,
                                 ref_id1 = ref_id1, ref_id2 = ref_id2)
  
  assign("ranked_genes", ranked_genes, .GlobalEnv)
  assign("gene_set", gene_set, .GlobalEnv)
  
  res <- fgsea(pathways = gene_set, stats = ranked_genes, 
               nperm = nperm, minSize = minSize, maxSize = maxSize)
  
  assign("gsea_res", res, .GlobalEnv)
  
  message("GSEA results for the comparison with the following parameters are stored in 'gsea_results'")
  message("Sample (left side of enrichment plot)")
  print(paste("samp_id1:", samp_id1, sep = " "))
  print(paste("samp_id2:", samp_id2, sep = " "))
  print(paste("marker1:", marker1, sep = " "))
  print(paste("marker2:", marker2, sep = " "))
  print(paste("gene_set:", deparse(substitute(gene_set)), sep = " "))
  
  message("Reference (right side of enrichment plot)")
  print(paste("ref_id1:", ref_id1, sep = " "))
  print(paste("ref_id2:", ref_id2, sep = " "))
  print(paste("marker1:", marker1, sep = " "))
  print(paste("marker2:", marker2, sep = " "))
  print(paste("gene_set:", deparse(substitute(gene_set)), sep = " "))
  message(cat("\n\n"))
  
  if (is.null(plot_individual)) {
    
    top_plotter(gsea_results = res, ranked_genes = ranked_genes, gene_set = gene_set,
                top_n = top_n, gseaParam = gseaParam)
  } else {
    
    hits <- c(grep(plot_individual, res$pathway, ignore.case = T, value = T))
    assign("hits", hits, .GlobalEnv)
    
    
    if (length(hits) > 1){
      multiple_hits <- t(t(hits))
      colnames(multiple_hits) <- "Multiple pathway matches"
      rownames(multiple_hits) <- c(1:length(hits))
      print(multiple_hits)
      n <- as.numeric(readline(prompt = "Multiple pathways are found. Select a number from the list above "))
      
      
      while(!n %in% 1:length(hits)) {
        n <- as.numeric(readline(prompt=paste0("Please pick a number between 1 and ", length(hits),":    ")))
      }
      assign("n", n , .GlobalEnv)
      
      annot_padj <- signif(as.numeric(res[res$pathway==hits[n], "padj"]), digits = 2)
      annot_NES <- signif(as.numeric(res[res$pathway==hits[n], "NES"]),digits=2)
      
      grob<- grobTree(textGrob(paste("adj.p: ", annot_padj, "\nNES:", annot_NES), x= 0.1, y=0.35, hjust = 0,
                               gp = gpar(col="red", fontsize=13, fontface="italic")))
      
      if(append_title ==F){
        
        plotEnrichment(pathway = gene_set[[hits[n]]], stats = ranked_genes) +
          labs(title = hits[n]) +
          annotation_custom(grob)
        
      } else {
        
        plotEnrichment(pathway = gene_set[[hits[n]]], stats = ranked_genes) +
          labs(title = paste0(hits[n], " (", samp_id1, "_&_", samp_id2, " vs ", ref_id1, "_&_", ref_id2,")")) +
          annotation_custom(grob)
        
      }
      

      
      ########################## 
      
    } else {
      
      annot_padj <- signif(as.numeric(res[res$pathway==hits, "padj"]), digits = 2)
      annot_NES <- signif(as.numeric(res[res$pathway==hits, "NES"]),digits=2)
      
      grob<- grobTree(textGrob(paste("adj.p: ", annot_padj, "\nNES:", annot_NES), x= 0.1, y=0.35, hjust = 0,
                               gp = gpar(col="red", fontsize=13, fontface="italic")))
      
      if(append_title ==F){
        
        plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
          labs(title = hits) +
          annotation_custom(grob)
        
      } else {
        
        plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
          labs(title = paste0(hits, " (", samp_id1, "_&_", samp_id2, " vs ", ref_id1, "_&_", ref_id2,")")) +
          annotation_custom(grob)
        
      }
      
    }
    
  }
  
}

##############################################################################################################
# # Examples of using master_gsea function
# master_gsea(exprs = exprs, marker1 = "ADGRE1", samp_id1 = "d12wt", ref_id1 = "d12ko",
#             gene_set = msig_GO ,top_n = 10, gseaParam = 0.3, plot_individual = "go_wound_healing")
# ggsave("d12_f480+_WTvsKO_woundhealing.png", width = 3.5, height = 2, units = "in", scale = 1.5)
# 
# 
# master_gsea(exprs = exprs, marker1 = "ADGRE1", samp_id1 = "d12wt", ref_id1 = "d12ko",
#             gene_set = msig_GO ,top_n = 10, gseaParam = 0.3, plot_individual = "gamma")
# ggsave("d12_f480+_WTvsKO_IFNg_response.png", width = 3.5, height = 2, units = "in", scale = 1.5)
# 
# 
# master_gsea(exprs = exprs, marker1 = "CD8A", marker2 = "CD3E", samp_id1 = "d12wt", ref_id1 = "d12ko",
#             gene_set = msig_GO ,top_n = 10, gseaParam = 0.3, plot_individual = "mitotic")
# ggsave("d12_cd8a+cd3e+_WTvsKO_pos_reg_mitosis.png", width = 3, height = 2, units = "in", scale = 1.5)
# 
# master_gsea(exprs = exprs, gene_set = msig_ALL, 
#             samp_id1 = "d12wt", samp_id2 = "Activated", 
#             ref_id1 = "d12wt", ref_id2 = "Act._Cycling", top_n = 10, gseaParam = 0.3)
#
#
# # After executing the code, you can use the following code to:
# # 1) Quickly examine enrichment results in table format (filtered based on p value and/or gene set of interest
# # 2) Visualize individual plots (for exploration)
# View(gsea_res[padj < 0.05, ])
# View(gsea_res[pathway %in% grep("neut", gsea_res$pathway, ignore.case = T, value = T) & padj <0.05,])
# plotEnrichment(pathway=msig_ALL[["GSE40666_NAIVE_VS_EFFECTOR_CD8_TCELL_DN"]], stats = ranked_genes)

##############################################################################################################
##############################################################################################################





