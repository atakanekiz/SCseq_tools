homolog_converter <- function(input_df, # data frame containing genes in rows, samples in columns. Rows and columns should be named
                              input_species = "mouse",
                              output_species = "human"){
  
  require(biomaRt)
  require(tibble)
  require(data.table)
  require(dplyr)
  
  # input_df <- dat
  
  input_df <- data.table(add_column(input_df, gene_name = rownames(input_df), .after = 0))
  
  genes_to_convert <- input_df$gene_name
  
  if(input_species == "mouse") inputmart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  if(output_species == "human") outputmart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  
  conversion_df = getLDS(filters = "mgi_symbol", 
                         values = genes_to_convert, 
                         attributes = c("mgi_symbol"), 
                         mart = inputmart,
                         attributesL = c("hgnc_symbol"), 
                         martL = outputmart, uniqueRows=T)
  
  input_df <- merge(conversion_df, input_df, by.x = "MGI.symbol", by.y="gene_name", sort = F)
  
  input_df <- input_df[,colnames(input_df) != "MGI.symbol"]
  
  input_df <- setDT(input_df)[,lapply(.SD, mean), by=HGNC.symbol] # MUCH FASTER than aggregate and dplyr functions
  
  human_genes <- input_df %>% pull(HGNC.symbol)  # Capture gene names to give it back as row names
  
  input_df <- as.data.frame(input_df)
  
  input_df <- input_df[, colnames(input_df) != "HGNC.symbol"] 
  
  rownames(input_df) <- human_genes
  
  input_df
  
}
