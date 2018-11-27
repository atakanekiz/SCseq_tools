# SCseq_tools
Functions to facilitate SCseq RNA seq data analysis


This repository features functions I commonly use in SCseq analysis.

The expanding list of functions:

- df_extractor
  - Tidy data frame that includes meta data (such as sample, cluster) and gene expression data
  - Rows constitute individual cells, columns are genes
  - Ability to converting mouse genes to human homologs

- gene_ranker
  - Allows user to select one of multiple metrics for ranking
  - Used in the gsea_plotter function for work under the hood

- top_plotter
  - Plots user-defined number of top GSEA results
  - Used in the gsea_plotter function for work under the hood

- gsea_plotter
  - Allows automatic subsetting of the data frame based on sample and cluster identity
  - User-defined number of top hits or a selected pathway can be shown
  - Regex-based name matching of pathway of interest that will prompt further user dialog in case of multiple hits

- Select gene plotter
  - Allows user to chose different plot types (Box, Column, Violin)
  - Uses tidy expression data frame mentioned above
  - Also shows pairwise statistics from all or pre-selected sample pairs
  - Allows subsetting the dataset based on sample and cluster identity
