gsea_plotter <- function(exprs = NULL, # Expression data frame (rows are cells, columns are genes and metadata)
                        pos_marker = NULL, # Genes to positively gate cells (cells expressing these markers will be subsetted)
                        neg_marker = NULL, # Genes to negatively gate cells (cells expressing these markers will be discarded)
                        sample_id, # Which cells will be selected as 'sample' (ie, the direction of rank ordering. Regex based string recognition
                        sample_cluster = NULL, # Cell clusters to include in analysis for the sample subset
                        reference_id, # Which cells will be selected as 'reference' (ie, the direction of rank ordering)
                        reference_cluster = NULL, # Cell clusters to include in analysis for the reference subset
                        method = "s2n", #"ttest", "difference", "ratio", "welch", "mwt", "bws"), # Method for ranking the genes. For MWT, see PMID: 18344518
                        humanize_exprs_df = T,
                        gene_set = "hallmark", #"go", "curated", "immune", "motif", "all", "custom"), 
                        nperm = 10000, 
                        minSize = 50, 
                        maxSize = 500,
                        top_n = 10, 
                        gseaParam = 1, 
                        plot_individual = NULL,
                        append_title = F){
  
  
  require(fgsea)
  require(dplyr)
  require(grid)
  require(ggplot2) 
  
  source("gene_ranker_v2.R")
  source("top_plotter.R")
  
  # Read molecular signatures database (MSigDB) gene lists. files must be stored 
  
  if(gene_set == "hallmark"){
    
    gene_set <- gmtPathways("msigdb_genelists/h.all.v6.2.symbols.gmt")
  
  } else if(gene_set == "go"){
    
    gene_set <- gmtPathways("msigdb_genelists/c5.all.v6.2.symbols.gmt")
      
  } else if(gene_set == "curated"){
    
    gene_set <- gmtPathways("msigdb_genelists/c2.all.v6.2.symbols.gmt")
      
  } else if(gene_set == "immune"){
    
    gene_set <- gmtPathways("msigdb_genelists/c7.all.v6.2.symbols.gmt")
      
  } else if(gene_set == "motif"){
    
    gene_set <- gmtPathways("msigdb_genelists/c3.all.v6.2.symbols.gmt")
      
  } else if(gene_set == "all"){
    
    gene_set <- gmtPathways("msigdb_genelists/msigdb.v6.2.symbols.gmt")
      
  } else {gene_set = gene_set} # can pass a named list composing of genes as character vectors in each list element

  if(is.null(reference_cluster)) {reference_cluster <- sample_cluster} # Specify only sample_cluster for subsetting on the same clusters in sample and reference
  
  
  ranked_genes <- gene_ranker(exprs = exprs,
                              pos_marker = pos_marker,
                              neg_marker = neg_marker,
                              sample_id = sample_id, 
                              sample_cluster = sample_cluster,
                              reference_id = reference_id,
                              reference_cluster = reference_cluster, 
                              method = method)
  
  # assign("ranked_genes", ranked_genes, .GlobalEnv)
  # assign("gene_set", gene_set, .GlobalEnv)
  
  

  
  res <- fgsea(pathways = gene_set, stats = ranked_genes, 
               nperm = nperm, minSize = minSize, maxSize = maxSize)
  
  # assign("gsea_res", res, .GlobalEnv)

  if (is.null(plot_individual)) {
    
    top_plotter(gsea_results = res, ranked_genes = ranked_genes, gene_set = gene_set,
                top_n = top_n, gseaParam = gseaParam)
  } else {
    
    hits <- c(grep(plot_individual, res$pathway, ignore.case = T, value = T))
    
    # assign("hits", hits, .GlobalEnv)
    
    
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
        
        arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = pos_marker, neg = neg_marker)
        select_non_null <- !sapply(arg_list, is.null)
        
          plot_subtitle <- paste(arg_list[select_non_null], names(arg_list[select_non_null]), sep="", collapse = "__")

        
        
        plotEnrichment(pathway = gene_set[[hits[n]]], stats = ranked_genes) +
          labs(title = paste0(hits[n], sample_id, " vs ", reference_id),
               subtitle = plot_subtitle)+
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
        
        arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = pos_marker, neg = neg_marker)
        select_non_null <- !sapply(arg_list, is.null)
        
        plot_subtitle <- paste(arg_list[select_non_null], names(arg_list[select_non_null]), sep="", collapse = "__")
        
        plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
          labs(title = paste0(hits[n], sample_id, " vs ", reference_id),
               subtitle = plot_subtitle)+
          annotation_custom(grob)
        
      }
      
    }
    
  }
  
}

##############################################################################################################
# # Examples of using master_gsea function

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# exprs <- readRDS("aging_exprs.rds")

exprs <- readRDS("aging_exprs.rds")

gsea_plotter(exprs, 
             sample_id = "Young \\(WT\\)", 
             reference_id = "Aged \\(WT\\)",
             sample_cluster = "NK",
             reference_cluster = "NK", gene_set = "hallmark")

















