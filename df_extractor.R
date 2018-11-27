# Make a function to create gene expression data frame to be used in GSEA calculations and gene plotting
###################################################################################################################################################################

df_extractor <- function(seurat_obj, 
                         metadata_to_extract = c(Cluster="ident",Sample="cond"), # User defined metadata column name that has sample information 
                         use_raw = F, # Set this to TRUE if you want to extract raw data
                         humanize = T # Convert mouse genes to human homologs?
){
  
  
  require(Seurat)
  require(tibble)
  require(dplyr)
  
  
  if(use_raw == F){
    
    exprs <- as.data.frame(as.matrix(seurat_obj@data))
    
  } else {
    
    exprs <- as.data.frame(as.matrix(seurat_obj@raw.data))
    }
  
  
  
  
  if(humanize==T){
    
    exprs <- add_column(exprs, gene = rownames(exprs), .after = 0)
    
    require(biomaRt)
    require(data.table)

    
    mouse_genes <- exprs$gene
    
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    convert_genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_genes , mart = mouse,
                           attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    
    exprs <- merge(convert_genes, exprs, by.x = "MGI.symbol", by.y="gene", sort = F)
    
    df <- exprs[,colnames(exprs) != "MGI.symbol"]
    
    df <- setDT(df)[,lapply(.SD, mean), by=HGNC.symbol] # MUCH FASTER than aggregate and dplyr functions
    
    human_genes <- df %>% pull(HGNC.symbol)  # Capture gene names to give it back as row names
    
    df <- as.data.frame(df)
    
    df <- df[, colnames(df) != "HGNC.symbol"] 
    
    rownames(df) <- human_genes
    
    
  } else {
    
    df <- exprs
    
  }
  
  
  # Transpose df
  
  df <- as.data.frame(t(df))
 
  for(i in names(metadata_to_extract)){
    
    if(metadata_to_extract[i] == "ident"){
      df <- add_column(df, !!i := slot(seurat_obj, metadata_to_extract[[i]]), .after = 0)
      next}
    
    df <- add_column(df, !!i := seurat_obj@meta.data[[metadata_to_extract[[i]]]], .after = 0)
    
    df <- add_column(df, Cell_id = colnames(seurat_obj@data))
    
  }
  
  df
  
}


################################################################################################

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# seurat_obj <- readRDS("c:/OConnell Lab/Experiments/June 2018/mir146ko-mir155tcko aging single cell sequencing/Seurat analyses/reanalysis_v2/combined_w_tsne.rds")



# # How to use:
# exprs <- df_extractor(seurat_obj = seurat_obj)


# saveRDS(exprs, "aging_exprs.rds")



