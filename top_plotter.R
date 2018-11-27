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