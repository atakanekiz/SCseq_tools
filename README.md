# SCseq_tools
Functions to facilitate SCseq RNA seq data analysis


This repository features functions I commonly use in SCseq analysis.

The capabilities include:

- Expression data frame extractor
  - Tidy data frame that includes meta data (such as sample, cluster) and gene expression data
  - Rows constitute individual cells, columns are genes
  - Ability to converting mouse genes to human homologs
  
- GSEA plotter
  - Gene ranker subfunction (multiple metrics)
  - Top GSEA results plotter subfunction (user-defined top hits)
  - Allows automatic subsetting of the data frame based on sample and cluster identity
  
- Select gene plotter
  - Allows user to chose different plot types (Box, Column, Violin)
  - Uses tidy expression data frame mentioned above
  - Also shows pairwise statistics from all or pre-selected sample pairs
  - Allows subsetting the dataset based on sample and cluster identity
