gsea_plotter <- function(exprs = NULL, # Expression data frame (rows are cells, columns are genes and metadata)
                        pos_marker = NULL, # Genes to positively gate cells (cells expressing these markers will be subsetted)
                        neg_marker = NULL, # Genes to negatively gate cells (cells expressing these markers will be discarded)
                        sample_id, # Which cells will be selected as 'sample' (ie, the direction of rank ordering. Regex based string recognition. Make sure you escape special characters such as parantheses
                        sample_cluster = NULL, # Cell clusters to include in analysis for the sample subset
                        reference_id, # Which cells will be selected as 'reference' (ie, the direction of rank ordering)
                        reference_cluster = NULL, # Cell clusters to include in analysis for the reference subset
                        method = "s2n", #"ttest", "difference", "ratio", "welch", "mwt", "bws"), # Method for ranking the genes. For MWT, see PMID: 18344518
                        gene_set = "hallmark", #"go", "curated", "immune", "motif", "all", "custom"), 
                        nperm = 10000, # number of permutations performed for calculations
                        minSize = 50, # minimum number of genes in a list to be considered in enrichment
                        maxSize = 500, # maximum number of genes in a list to be considered in enrichment
                        top_n = 10, # Report top enrichment results
                        gseaParam = 1, # Changes the size of gene bars in summary plots
                        plot_individual = NULL, #
                        append_title = F, # add groups in comparison to as subtitle in individual plots
                        top_plots_title = T, # Plotting title to summarize the groups in comparison
                        seed = 123, # random seed
                        keep_results = T, # Save enrichment results as a global object?
                        save_png = F, 
                        png_units = "in", 
                        png_width = 4, # good size for individual plots. For summary table increase the size
                        png_height = 3,
                        append_to_filename = "" # add a custom string to the png filename
                        ){
  
  set.seed(seed)
  
  require(fgsea)
  require(dplyr)
  require(grid)
  require(ggplot2) 
  
  # source("gene_ranker_v2.R")
  # source("top_plotter.R")
  
  # Read molecular signatures database (MSigDB) gene lists. files must be stored 
  
  if(gene_set == "hallmark"){
    
    gene_set <- hallmark
  
  } else if(gene_set == "go"){
    
    gene_set <- go
      
  } else if(gene_set == "curated"){
    
    gene_set <- curated
      
  } else if(gene_set == "immune"){
    
    gene_set <- immune
      
  } else if(gene_set == "motif"){
    
    gene_set <- motif
      
  } else if(gene_set == "all"){
    
    gene_set <- all
      
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
  
  if(keep_results) assign("ranked_genes", ranked_genes, .GlobalEnv)
  # assign("gene_set", gene_set, .GlobalEnv)
  
  

  
  res <- fgsea(pathways = gene_set, stats = ranked_genes, 
               nperm = nperm, minSize = minSize, maxSize = maxSize)
  
  if(keep_results) assign("gsea_res", res, .GlobalEnv)

  if (is.null(plot_individual)) {
    
    if(top_plots_title == T) {
      
      arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker, collapse = "."))
      select_non_null <- !sapply(arg_list, function(x) {identical(x, "")})
      
      
     
      main_title <- paste(sample_id, "vs", reference_id)
      plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null],  sep=": ", collapse = "__")
      
      plot_title <- paste0(main_title,"\n", plot_subtitle)
      
    } else { plot_title = ""}
    
    plot_grob <- top_plotter(gsea_results = res, ranked_genes = ranked_genes, gene_set = gene_set,
                top_n = top_n, gseaParam = gseaParam, plot_title=plot_title, do.plot=F)
    grid.arrange(plot_grob)
    
    
  } else {
    
    hits <- c(grep(plot_individual, res$pathway, ignore.case = T, value = T))
    
    # assign("hits", hits, .GlobalEnv)
    
    
    if (length(hits) > 1){
      multiple_hits <- t(t(hits))
      colnames(multiple_hits) <- "Multiple pathway matches"
      rownames(multiple_hits) <- c(1:length(hits))
      print(multiple_hits)
      num <- as.numeric(readline(prompt = "Multiple pathways are found. Select a number from the list above "))
      
      
      while(!num %in% 1:length(hits)) {
        num <- as.numeric(readline(prompt=paste0("Please pick a number between 1 and ", length(hits),":    ")))
      }
      assign("num", num , .GlobalEnv)
      
      annot_padj <- signif(as.numeric(res[res$pathway==hits[num], "padj"]), digits = 2)
      annot_NES <- signif(as.numeric(res[res$pathway==hits[num], "NES"]),digits=2)
      annot_ES <- signif(as.numeric(res[res$pathway==hits[num], "ES"]),digits=2)
      
      # grob<- grobTree(textGrob(paste("adj.p: ", annot_padj, "\nNES:", annot_NES), 
      #                          x= 0.1, y=annot_ES, 
      #                          hjust = 0,
      #                          gp = gpar(col="red", fontsize=3, fontface="italic")))
      
      annot_text <- paste("adj.p: ", annot_padj, "\nNES:", annot_NES)
      
      if(append_title ==F){
        
        # plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) +
        #   labs(title = hits[num]) +
        #   annotation_custom(grob)
        
        plot_grob <- plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) +
          labs(title = hits[num]) +
          annotate("text", x=x_pos, y=annot_ES/2, label = annot_text, colour = "red", size=4)+
          theme(plot.title = element_text(size=5, hjust = 0.5))
        print(plot_grob)
        
      } else {
        
        arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = pos_marker, neg = neg_marker)
        select_non_null <- !sapply(arg_list, function(x) {identical(x, "")})
        select_non_null2 <- !sapply(arg_list, is.null)
        select_non_null <- as.logical(select_non_null * select_non_null2)
        
        plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep=": ", collapse = "__")

        
        
        # plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) +
        #   labs(title = paste0(hits[num], sample_id, " vs ", reference_id),
        #        subtitle = plot_subtitle)+
        #   annotation_custom(grob)
        
        plot_grob <- plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) +
          labs(title = hits[num],
               subtitle = paste(sample_id, "vs", reference_id, plot_subtitle))+
          annotate("text" , x=x_pos, y=annot_ES/2, label = annot_text, colour = "red", size=4)+
          theme(plot.title = element_text(size=10, hjust = 0.5),
                plot.subtitle = element_text(size=6, hjust = 0.5))
        print(plot_grob)
        
      }
      
      
      
      ########################## 
      
    } else {
      
      annot_padj <- signif(as.numeric(res[res$pathway==hits, "padj"]), digits = 2)
      annot_NES <- signif(as.numeric(res[res$pathway==hits, "NES"]),digits=2)
      annot_ES <- signif(as.numeric(res[res$pathway==hits, "ES"]),digits=2)
      x_pos <- length(ranked_genes)/5
      
      # grob<- grobTree(textGrob(paste("adj.p: ", annot_padj, "\nNES:", annot_NES),
      #                          x= 0.1, y=annot_ES, 
      #                          hjust = 0,
      #                          gp = gpar(col="red", fontsize=3, fontface="italic")))
      
      annot_text <- paste("adj.p: ", annot_padj, "\nNES:", annot_NES)
      
      
      if(append_title ==F){
        
        # plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
        #   labs(title = hits) +
        #   annotation_custom(grob)
        
        plot_grob <- plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
          labs(title = hits) +
          annotate("text", x=x_pos, y=annot_ES/2, label = annot_text, colour = "red", size=4)
        print(plot_grob)
        
      } else {
        
        arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker, collapse = "."))
        
        select_non_null <- !sapply(arg_list, function(x) {identical(x, "")})
        select_non_null2 <- !sapply(arg_list, is.null)
        select_non_null <- as.logical(select_non_null * select_non_null2)
        
        plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null],  sep=": ", collapse = " ")
        
        # plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
        #   labs(title = paste0(hits, sample_id, " vs ", reference_id),
        #        subtitle = plot_subtitle)+
        #   annotation_custom(grob)
        
        plot_grob <- plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) +
          labs(title = hits,
               subtitle = paste(sample_id, "vs", reference_id, plot_subtitle))+
          annotate("text", x=x_pos, y=annot_ES/2, label = annot_text, colour = "red", size=4)+
          theme(plot.title = element_text(size=10, hjust = 0.5),
                plot.subtitle = element_text(size=6, hjust = 0.5))
        print(plot_grob)
      }
      
    }
    
  }
  
  if(save_png == T){
    
   sample_cluster <- str_replace_all(sample_cluster, "[:punct:]|[:space:]", "")
   reference_cluster <-  str_replace_all(reference_cluster, "[:punct:]|[:space:]", "")
   sample_id  <-  str_replace_all(sample_id, "[:punct:]|[:space:]", "")
   reference_id <-   str_replace_all(reference_id, "[:punct:]|[:space:]", "")
    
    arg_list <- list(SAMPclus = sample_cluster, 
                     REFclus = reference_cluster, 
                     SAMPid = sample_id, 
                     REFid = reference_id,
                     POSmarker = paste(pos_marker, collapse = "."), 
                     NEGmarker = paste(neg_marker, collapse = "."))
    
    select_non_null <- !sapply(arg_list, function(x) {identical(x, "")})
    select_non_null2 <- !sapply(arg_list, is.null)
    select_non_null <- as.logical(select_non_null * select_non_null2)
    

    
    
    if(sum(select_non_null) == 0) {
      
      filename <- paste0("Unsubsetted_data_plots__", ".pdf")
      
    } else {
      
      filename <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep="-", collapse = "___")
      filename <- paste0(filename, "___", append_to_filename, ".png")
    }
    
    ggsave(plot = plot_grob, filename = filename, width = png_width, height = png_height, units = png_units)
  
  } #else grid.arrange(plot_grob)
  
}

##############################################################################################################
# # Examples of using master_gsea function

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# exprs <- readRDS("aging_exprs.rds")

# exprs <- readRDS("aging_exprs.rds")

# gsea_plotter(exprs, 
#             sample_id = "Young \\(WT\\)", 
#             reference_id = "Aged \\(WT\\)",
#             sample_cluster = "NK",
#             reference_cluster = "NK", gene_set = "hallmark")

















