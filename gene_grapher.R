# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# exprs <- readRDS("aging_exprs.rds")


gene_grapher <- function(exprs, # expression dataframe (generated with df_extractor function)
                         genes_to_plot, # genes to plot 
                         x_variable = c("Sample", "Cluster"), # what to show on x-axis
                         clusters_to_plot = NULL, # select clusters to be subsetted for plot
                         pos_marker = NULL, # select cells expressing these markers. Can provide a character vector of length 1 or more
                         neg_marker = NULL, # select cells not expressing these markers. Can provide a character vector of length 1 or more
                         plot_type = c("box", "bar", "violin"), # the type of plot to be generated
                         add_jitter = T, # Add transparent jitter data points
                         add_mean = T, # Add a red colored point indicating mean value
                         add_median = T, # Add a blue colors point indicating median value
                         sort_plots = F, # Alphabetical ordering of plots based on gene name
                         colors_to_use = NULL, # Default is rainbow palette. You can provide a character vector
                         show_stats = T, # Calculate and show statistics on graph?
                         comparisons = NULL, # Which stats to show. A list of character vectors (pair-wise)
                         stat_method = c("t.test", "wilcox.test"), # Use parametric t-test or nonparametric wilcoxon test
                         p_adj_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"), # How to adjust p-values
                         save_pdf = T, # Save results in a pdf file in the workspace.
                         append_to_filename = "",
                         output_plot = F, # Set to TRUE if you'd like to show graph in RMD or viewer
                         assign_global_plotlist = F, # Stores plot_list as a global object for further manipulations
                         show_progress = F, # Returns a message to the console to let you know which gene is being plotted 
                         image_columns = 4, # How many columns of images should be on one page
                         image_rows = 4, # How many rows of images should be on one page
                         ... # Not used currently, but can be included to pass other arguments to graphing functions
)


  {
  
  require(ggpubr)
  require(dplyr)
  require(stringr)
  require(purrr)
  require(colorspace)
  require(gtools)
  
  colnames(exprs) <- str_replace_all(colnames(exprs), "[:punct:]|[:space:]", "_")
  genes_to_plot <- str_replace_all(genes_to_plot, "[:punct:]|[:space:]", "_")
  genes_to_plot <- genes_to_plot[genes_to_plot %in% colnames(exprs)]
  

  
  if(sort_plots ==T) {genes_to_plot <- mixedsort(genes_to_plot)}
  
  # Subset clusters of interest
  if(!is.null(clusters_to_plot)) {
    clusters_to_plot <- gsub("\\+", "\\\\+", clusters_to_plot) # Escape special characters for correct subsetting
    clusters_to_plot <- gsub("\\-", "\\\\-", clusters_to_plot)
    exprs <- filter(exprs, str_detect(Cluster, regex(clusters_to_plot, ignore_case=T, comments = F)))
    }
  
  # Subset cells that express markers of interest at any level (>0)
  if(!is.null(pos_marker)) {
    
    # pos_marker_split <- trimws(strsplit(pos_marker, split = ",")[[1]])
    # exprs <- filter_(exprs, paste(pos_marker_split , "!= 0", collapse = " & "))   # Can be used to pass a comma separated long string "gene1 , gene2"
    
    exprs <- filter_(exprs, paste(pos_marker , "!= 0", collapse = " & "))
    
    }
  
  
  # Discard cells that express a marker gene which we want to negatively gate in our analyses
  if(!is.null(neg_marker)) {
    
    # neg_marker_split <- trimws(strsplit(neg_marker, split = ",")[[1]])
    # exprs <- filter_(exprs, paste(neg_marker_split , "== 0", collapse = " & "))  # Can be used to pass a comma separated long string "gene1 , gene2"
    
    exprs <- filter_(exprs, paste(neg_marker_split , "== 0", collapse = " & "))
    
    }
  
  
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
  
  plot_list <- list() # Initialize plot_list to store plots to be printed
  
  for(gene in genes_to_plot){
    
  if(show_progress == T) {message(paste("Plotting", gene))}  

  if(show_stats == T){
    
    tryCatch(error = function(x){warning(paste("Statistics cannot be computed for ", gene))},
             {
    
    # Compute p-values
    comparison_formula <- paste0(gene, "~", x_variable) %>%
      as.formula()
    stat_test <- compare_means(
      comparison_formula,  data = exprs,
      method = stat_method,
      p.adjust.method = p_adj_method
    )
    
    # If a comparison list is provided, extract the comparisons of interest for plotting
    if(!is.null(comparisons)){
      
      stat_test <- do.call(rbind, comparisons) %>% # rbind to a matrix
        as.data.frame %>% # convert to a data.frame
        set_names(c("group1", "group2")) %>% # change the column names
        inner_join(stat_test) # and inner join
      
    }
    
    # P-value y coordinates
    y_max <- exprs %>% 
      pull(gene) %>% max(na.rm = TRUE) 
    
    
    p_value_y_coord <- rep(y_max*5/6, nrow(stat_test))
    
    step_increase <- (1:nrow(stat_test))*(y_max/5)
    p_value_y_coord <- p_value_y_coord + step_increase
    
    stat_test <- stat_test %>%
      mutate(
        y.position =  p_value_y_coord, 
        p.adj = format.pval(p.adj, digits = 1)
      )
    
}) # close tryCatch
    
  } # close if(show_stats=T)
    
    if(is.null(colors_to_use)){
      
      group_number <- length(levels(as.factor(pull(exprs, x_variable))))
      colors_to_use <- rainbow_hcl(group_number)
      
    }
    
    if(grepl("^[0-9]", gene)) gene <- paste0("`", gene, "`")  # Append backticks for colnames starting with number
    
    if(plot_type == "box"){
      
      
      
      p <- ggboxplot(exprs, x= x_variable, y = gene, fill = x_variable, palette = colors_to_use, 
                     title = gene, ylab="Expression", outlier.shape=NA, xlab = F)+
        rremove("legend")+
        rotate_x_text(angle=45)
      
      } else if(plot_type == "bar"){
        
        
      
      p <- ggbarplot(exprs, x= x_variable, y = gene, fill = x_variable, palette = colors_to_use,
                     title = gene, ylab="Expression",
                     add =  "mean_se" , xlab = F)+
        rremove("legend")+
        rotate_x_text(angle=45)
      
     
      
    } else if(plot_type == "violin"){
      
      
      
      p <- ggviolin(exprs, x= x_variable, y = gene, fill = x_variable, palette = colors_to_use,
               title = gene, ylab="Expression",  scale="width", trim = T)+ #draw_quantiles = c(0.25, 0.5, 0.75) can be added for marking IQR and median
        rremove("legend")+
        rotate_x_text(angle=45)
      
      
    }
    
    if(add_jitter == T){ p <- ggadd(p, add = "jitter", alpha=0.2, size = 0.2) }
    
    if(add_mean == T){ p <- ggadd(p, add = "mean", color = "red", size = 0.3)}
    
    if(add_median == T){ p <- ggadd(p, add = "median", color = "blue", size = 0.3)}
    
    if(show_stats == T){ tryCatch(error=function(x){},
                                  
                                  {
      p <- p + expand_limits(y=y_max*2) +
        stat_pvalue_manual(stat_test, label = "p.signif", size = 2.5)
      
                                  })
    }
    
   plot_list[[gene]] <- p 
    
    
  }
  
  if(assign_global_plotlist == T) { assign("plot_list", plot_list, envir = .GlobalEnv)}
    
  if(output_plot ==T) print(ggarrange(plotlist = plot_list, ncol = image_columns, nrow = image_rows))
  
  if(save_pdf == T){
    
    arg_list <- list(cluster = clusters_to_plot, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker, collapse = "."))
    select_non_null <- !sapply(arg_list, is.null)
    
    if(sum(select_non_null) == 0) {
      
      filename <- paste0("Unsubsetted_data_plots__", stat_method, ".pdf")
        
    } else {
      
      filename <- paste(arg_list[select_non_null], names(arg_list[select_non_null]), collapse = "__")
      filename <- paste0(filename,"__", stat_method, append_to_filename, ".pdf")
    }
    

      
    ggexport(plotlist = plot_list, filename = filename,
             nrow = image_rows, ncol=image_columns, 
             height = image_rows*2.5, width = image_columns*1.8, res = 300)
    
  } 
    
  }
  
  

##############################################################################################
###################   Test   #################################################################
##############################################################################################


# ms <- read.delim("mouse_immune_process.txt", header = F)
# 
# plot_clust <- c("Nai Cd8+", "Nai Cd4+") #"Act Cd8+", "Mem Cd4+"
# 
# 
# for(i in plot_clust) {
# 
# 
# 
# 
# 
# 
# gene_grapher(exprs = exprs,
#              genes_to_plot = ms$V1,
#              x_variable = "Sample",
#              clusters_to_plot = "Mem Cd4+",
#              pos_marker = NULL,
#              neg_marker = NULL,
#              plot_type = "box",
#              add_jitter = T,
#              add_mean = T,
#              add_median = F,
#              colors_to_use = c("dodgerblue", "coral2", "gray60", "gold2"),
#              show_stats = T,
#              comparisons = list(c("Young (WT)", "Aged (WT)"), 
#                                 c("Aged (WT)", "Aged (146 KO)"), 
#                                 c("Aged (146 KO)", "Aged (DKO)"), 
#                                 c("Aged (WT)", "Aged (DKO)")),
#              stat_method = "wilcox.test",
#              p_adj_method = "holm",
#              save_pdf = T,
#              append_to_filename = "_immune_process_GO",
#              output_plot = F,
#              assign_global_plotlist = F,
#              show_progress = T,
#              image_columns = 4,
#              image_rows = 4)
# 
#              
#  }

