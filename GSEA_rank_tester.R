set.seed(10)

gene_num <- 100

test <- data.frame(sample1=c(rnorm(gene_num/2, mean=50, sd=10), rnorm(gene_num/2, mean=20, sd=10)),
                   sample2=c(rnorm(gene_num/2, mean=50, sd=10), rnorm(gene_num/2, mean=20, sd=10)),
                   sample3=c(rnorm(gene_num/2, mean=50, sd=10), rnorm(gene_num/2, mean=20, sd=10)),
                   sample_a4=c(rnorm(gene_num/2, mean=50, sd=10), rnorm(gene_num/2, mean=20, sd=10)),
                   reference_b1=c(rnorm(gene_num/2, mean=20, sd=10), rnorm(gene_num/2, mean=40, sd=10)),
                   reference_b2=c(rnorm(gene_num/2, mean=20, sd=10), rnorm(gene_num/2, mean=40, sd=10)),
                   reference_b3=c(rnorm(gene_num/2, mean=20, sd=10), rnorm(gene_num/2, mean=40, sd=10)))



genes_row <- character()

for(i in 1:gene_num){
  
  add_gene <- paste(sample(LETTERS,4), collapse = "")
  
  genes_row <- c(genes_row, add_gene)
  
  }

rownames(test) <- genes_row


pheatmap::pheatmap(test)




# S2N ranking
##############################################################################
grp1 <- test[,1:4]
grp2 <- test[,5:7]


grp1_means <- rowMeans(grp1)
grp2_means <- rowMeans(grp2)

grp1_sd <- apply(grp1, 1, sd)
grp2_sd <- apply(grp2, 1, sd)

s2n_ranked <- (grp1_means - grp2_means) / (grp1_sd + grp2_sd)

s2n_ranked <- sort(s2n_ranked, decreasing = T)
# s2n_ranked



# tTest ranking
##############################################################################
grp1 <- test[,1:4]
grp2 <- test[,5:7]


grp1_means <- rowMeans(grp1)
grp2_means <- rowMeans(grp2)

grp1_sd <- apply(grp1, 1, sd)
grp2_sd <- apply(grp2, 1, sd)

grp1_n <- dim(grp1)[2]
grp2_n <- dim(grp2)[2]

ttest_ranked <- (grp1_means - grp2_means) / sqrt((grp1_sd^2/grp1_n)+(grp2_sd^2/grp2_n))
ttest_ranked <- sort(ttest_ranked, decreasing = T)
# ttest_ranked


# Ratio ranking (Recommended for log scale data)
##############################################################################
grp1 <- test[,1:4]
grp2 <- test[,5:7]


grp1_means <- rowMeans(grp1)
grp2_means <- rowMeans(grp2)

ratio_ranked <- grp1_means / grp2_means
ratio_ranked <- sort(ratio_ranked, decreasing = T)
# ratio_ranked


# Difference ranking (Recomended for natural scale data)
##############################################################################
grp1 <- test[,1:4]
grp2 <- test[,5:7]


grp1_means <- rowMeans(grp1)
grp2_means <- rowMeans(grp2)

difference_ranked <- grp1_means - grp2_means
difference_ranked <- sort(difference_ranked, decreasing = T)
# difference_ranked


# log2ratio ranking (Recomended for natural scale data)
##############################################################################
grp1 <- test[,1:4]
grp2 <- test[,5:7]


grp1_means <- rowMeans(grp1)
grp2_means <- rowMeans(grp2)

log2ratio_ranked <- log((grp1_means / grp2_means),2)
log2ratio_ranked <- sort(log2ratio_ranked, decreasing = T)
# log2ratio_ranked



##############################################################################
############### Use GeneSelector package for ranking #########################
##############################################################################

library(GeneSelector)

GeneSelector_tstat <- RankingTstat(x= as.matrix(test), y = factor(c(rep(1, times=4), rep(2, times=3))))
GeneSelector_tstat_res <- sort(GeneSelector_tstat@ranking, decreasing = F)

GeneSelector_welch <- RankingWelchT(x= as.matrix(test), y = factor(c(rep(1, times=4), rep(2, times=3))))
GeneSelector_welch_res <- sort(GeneSelector_welch@ranking, decreasing = F)

GeneSelector_difference <- RankingFC(x= as.matrix(test), y = factor(c(rep(1, times=4), rep(2, times=3))))
GeneSelector_difference_res <- sort(GeneSelector_difference@ranking, decreasing = F)


##############################################################################
############### Use mwt package for ranking ##################################
##############################################################################


library(mwt)

mwt_res <- mwt(as.matrix(test), grp = factor(c(1,1,1,1,2,2,2)))



##############################################################################
############### Use bws package for ranking ##################################
##############################################################################
library(BWStest)

bws_ranked <- numeric()

for(i in rownames(test)){
  
  samp <- as.numeric(test[i,1:4])
  ref <- as.numeric(test[i,5:7])
  
  bws_ranked[[i]] <- bws_stat(samp, ref)
  
}




##############################################################################
############################# Visualization ##################################
##############################################################################

# Compare results

res_list<-list(GeneSelector_difference = GeneSelector_difference@statistic, 
               GeneSelector_welch = GeneSelector_welch@statistic,
               GeneSelector_tstat = GeneSelector_tstat@statistic, 
               difference_ranked = difference_ranked, 
               ttest_ranked = ttest_ranked, 
               s2n_ranked= s2n_ranked,
               ratio_ranked = ratio_ranked, 
               log2ratio_ranked = log2ratio_ranked,
               mwt_res = mwt_res$MWT,
               bws_ranked = bws_ranked)

combined_res <- do.call(cbind.data.frame, res_list)


library(corrplot)
library(RColorBrewer)
M <-cor(combined_res)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=10, name="RdYlBu"), method = "pie")


# Show ties in BWS statistics

sum(duplicated(bws_ranked))

