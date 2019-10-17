# Suitable to extract data frame from Seurat v3 object


# Make a function to create gene expression data frame to be used in GSEA calculations and gene plotting
###################################################################################################################################################################

df_extractor <- function(seurat_obj, 
                         metadata_to_extract = c(Cluster="orig_clusters",Sample="group"), # User defined metadata column name that has sample information 
                         # use_raw = F, # Set this to TRUE if you want to extract raw data
                         humanize = T, # Convert mouse genes to human homologs?
                         assay = NULL, # which assay data to use. If null, default assay will be used
                         slot = NULL  # which specific information to pull. If null, data slot will be used

                         ){
  
  
  require(Seurat)
  require(tibble)
  require(dplyr)
  
  
 if(is.null(assay)) assay = DefaultAssay(seurat_obj)
 if(is.null(slot)) slot = "data"
    
    exprs <- as.data.frame(as.matrix(GetAssayData(object = seurat_obj, 
                                                  assay = assay,
                                                  slot = slot)))
    
    
  
  
  
  
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
  
  # df_org <- df
  # df <- df_org
  
  # Transpose df
  
  df <- as.data.frame(t(df))
 
  if(is.null(names(metadata_to_extract))) {
    
    warning("'metadata_to_extract' argument is an unnamed vector. Column names are created automatically. To ensure correct functioning of other downstream functions, ensure the column names are appropriate.")
    names(metadata_to_extract) <- metadata_to_extract
    
  }
  
  for(i in names(metadata_to_extract)){
    
      df <- add_column(df, !!i := seurat_obj[[metadata_to_extract[[i]]]][[1]], .after = 0)
  }
  
  df
  
}


################################################################################################

# # How to use:
# exprs <- df_extractor(seurat_obj = seurat_obj)



