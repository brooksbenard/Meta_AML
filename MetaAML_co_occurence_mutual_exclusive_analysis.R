# MetaAML_co_occurence_mutual_exlclusive_analysis
# Brooks Benard
# bbenard@stanford.edu
# 11/25/2019
# This script impliments several co-occuring and mutually exclusive mutation analysis packages on MetaAML data

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('tydyr')) install.packages('tydyr'); library('tydyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('cometExactTest')) install.packages('cometExactTest'); library('cometExactTest')
library(discover)
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('maditr')) install.packages('maditr'); library('maditr')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('epitools')) install.packages('epitools'); library('epitools')
library(corrplot)
library(RColorBrewer)

# create the data matrix to use in the function
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

final_data_matrix_2_sub <- subset(final_data_matrix, Subset == "de_novo")

# make sure that the FLT3 symbols are the same
for(i in 1:nrow(final_data_matrix_2_sub)){
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "ITD"){
    final_data_matrix_2_sub$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "SNV"){
    final_data_matrix_2_sub$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "Deletion"){
    final_data_matrix_2_sub$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "INDEL"){
    final_data_matrix_2_sub$Gene[i] <- "FLT3-ITD"
  }
}

# elect only recurrent mutations to perform the anlaysis on
final_data_matrix_2_sub$mut_freq_gene = NULL
pt_gene <- as.data.frame(unique(final_data_matrix_2_sub[,c(1:2)]))
pt_gene$mut_freq_gene <- NA

for(i in 1:nrow(pt_gene)){
  pt_gene_sub <- subset(pt_gene, pt_gene$Gene == pt_gene[i,2])
  pt_gene$mut_freq_gene[i] <- as.numeric(nrow(pt_gene_sub))
}
final_data_matrix_2_sub <- final_data_matrix_2_sub %>% left_join(pt_gene, by=c("Sample","Gene"))


temp1 <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$mut_freq_gene >= 25)
# temp1 = final_data_matrix_2_sub
# select informative columns
final_data_matrix_2_sub = select(final_data_matrix_2_sub, Sample, Gene, VAF)

# transform the data frame
final_data_matrix_2_sub <- dcast(final_data_matrix_2_sub, Sample ~ Gene, value.var="VAF")

rownames(final_data_matrix_2_sub) <- final_data_matrix_2_sub$Sample
final_data_matrix_2_sub$Sample <- NULL
final_data_matrix_2_sub[final_data_matrix_2_sub != 0] <- 1
temp_for_odds_ratio = final_data_matrix_2_sub
final_data_matrix_2_sub <- as.data.frame(t(final_data_matrix_2_sub))

# run the DISCOVER analysis
events <- discover.matrix(final_data_matrix_2_sub)
subset <- rowSums(final_data_matrix_2_sub) > 25
result.mutex <- pairwise.discover.test(events[subset, ])
result.mutex
print(result.mutex, fdr.threshold=0.05)
result.mutex = as.data.frame(result.mutex)




# odds ratio and fisher's exact p-value for each interaction ####
# find all pairwise interactions
genes=as.data.frame(t(combn(as.vector(unique(temp1$Gene)),2)))

n=1
results_list = list()

for(i in 1:nrow(genes)){
  print(i)
  gene_1=as.character(genes[i,1])
  gene_2=as.character(genes[i,2])
  
  temp = select(temp_for_odds_ratio, gene_1, gene_2)
   temp = as.data.frame(table(temp)) 
  
  # calculate odds ratio and fishers exact p-value for co-occurence
   # fishers exact test
   fish_table = data.frame(matrix(NA, nrow = 2, ncol = 2))
   
   fish_table[1,1] = temp$Freq[1]
   fish_table[1,2] = temp$Freq[2]
   fish_table[2,1] = temp$Freq[3]
   fish_table[2,2] = temp$Freq[4]
   
   fishers_p = fisher.test(fish_table)$p.value
   
   # odds ratio
   odds_table=as.matrix(fish_table)
   
   # oddsratio function breaks when there are no co-occurences so manually check that this isn't the case and correct
   odds_ratio = (temp$Freq[1]*temp$Freq[4])/(temp$Freq[2]*temp$Freq[3])
   if(odds_ratio != 0){
     odds_ratio = oddsratio(odds_table, conf.level = 0.95)$measure[2]
   } else if(odds_ratio == 0){
     # because plotting cases with an odds of zero really messes up the corrplot interpretability, set these to the next lowest odds ratio in the dataset
     odds_ratio = 0.04915185
   }

  # store the results in a dataframe
  odds <- data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(odds) <- c("gene1", "gene2", "odds_ratio", "fishers_exact")
  
  odds[1,1] <- gene_1
  odds[1,2] <- gene_2
  odds[1,3] <- odds_ratio
  odds[1,4] <- fishers_p
  
  results_list[[n]] <- odds
  n=n+1   
}
temp_final_1 = as.data.frame(do.call(rbind, results_list))

n=1
results_list = list()

for(i in 1:nrow(genes)){
  gene_1=as.character(genes[i,2])
  gene_2=as.character(genes[i,1])
  
  temp = select(temp_for_odds_ratio, gene_1, gene_2)
  temp = as.data.frame(table(temp))
  
  # calculate odds ratio and fishers exact p-value for co-occurence
# fisher's exact
  fish_table = data.frame(matrix(NA, nrow = 2, ncol = 2))
  
  fish_table[1,1] = temp$Freq[1]
  fish_table[1,2] = temp$Freq[2]
  fish_table[2,1] = temp$Freq[3]
  fish_table[2,2] = temp$Freq[4]
  
  fishers_p = fisher.test(fish_table)$p.value
  
  # odds ratio
  odds_table=as.matrix(fish_table)
  
  # oddsratio function breaks when there are no co-occurences so manually check that this isn't the case and correct
  odds_ratio = (temp$Freq[1]*temp$Freq[4])/(temp$Freq[2]*temp$Freq[3])
  if(odds_ratio != 0){
    odds_ratio = oddsratio(odds_table, conf.level = 0.95)$measure[2]
  } else if(odds_ratio == 0){
    # because plotting cases with an odds of zero really messes up the corrplot interpretability, set these to the next lowest odds ratio in the dataset
    odds_ratio = 0.04915185
  }

  
  # store the results in a dataframe
  odds <- data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(odds) <- c("gene1", "gene2", "odds_ratio", "fishers_exact")
  
  odds[1,1] <- gene_1
  odds[1,2] <- gene_2
  odds[1,3] <- odds_ratio
  odds[1,4] <- fishers_p
  
  results_list[[n]] <- odds
  n=n+1   
}
temp_final_2 = as.data.frame(do.call(rbind, results_list))

temp_final = unique(rbind(temp_final_1, temp_final_2))
# use temp_final later in visualizing the fishers exact results


# combine the DISCOVER results with the odds ratio numbers
DISCOVER_final = left_join(temp_final, result.mutex, by = c("gene1", "gene2"))

temp_final_odds = select(DISCOVER_final, gene1, gene2, odds_ratio)
temp_final_odds <- dcast(temp_final, gene1 ~ gene2, value.var="odds_ratio")
temp_final_odds = as.data.frame(temp_final_odds)
rownames(temp_final_odds) <- as.character(temp_final_odds$gene1)
temp_final_odds$gene1 <- NULL

# res1 <- cor.mtest(temp_final_odds, conf.level = .95)

for(i in 1:nrow(temp_final_odds)){
  for(j in 1:ncol(temp_final_odds)){
    # if( !is.na(temp_final_odds[i,j]) & temp_final_odds[i,j] == 0){
    #   temp_final_odds[i,j] = NA
    # }
    if(!is.na(temp_final_odds[i,j]) & temp_final_odds[i,j] != 0){
      temp_final_odds[i,j] = log(temp_final_odds[i,j])
    } 
  }
}

temp_final_odds = as.matrix(temp_final_odds)
temp_final_odds[is.na(temp_final_odds)] <- 0

# need to make sure the q-values calculated from the DISCOVER method are the same for each pait
result.mutex_2 = result.mutex[,c(2,1,3,4)]
colnames(result.mutex_2)[1:2] = c("gene1", "gene2")
DISCOVER_final_2 = left_join(DISCOVER_final, result.mutex_2, by = c("gene1", "gene2"))
DISCOVER_final_2$q.value = dplyr::coalesce(DISCOVER_final_2$q.value.x, DISCOVER_final_2$q.value.y)

# write out results file
write.csv(DISCOVER_final_2, "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/DISCOVER_results.csv", row.names=FALSE)



temp_final_q = select(DISCOVER_final_2, gene1, gene2, q.value)
temp_final_q <- reshape2::dcast(temp_final_q, gene1 ~ gene2, value.var="q.value")
temp_final_q = as.data.frame(temp_final_q)
rownames(temp_final_q) <- temp_final_q$gene1
temp_final_q$gene1 <- NULL
temp_final_q[is.na(temp_final_q)] <- 1
temp_final_q = as.matrix(temp_final_q)

pdf(file = "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/MetaAML_mutation_correlation_de_novo.pdf", width = 7.5, height = 7.5)

corrplot(temp_final_odds, is.corr = F, type="upper", order="hclust",tl.col="black", outline = F, addgrid.col = "lightgrey",
         col = brewer.pal(n = 8, name = "RdBu"), diag=FALSE, p.mat = temp_final_q, insig = "label_sig",
         sig.level = c(.001, .01, .1), pch.cex = .9, pch.col = "black", na.label = "square", na.label.col = "white")
dev.off()





# Fisher's exact test  ####
# plot results from Fisher's Exact analysis
# log transform the odds ratio so that it can be easily visualized in the corrplot
temp_final$odds_ratio_log = log(temp_final$odds_ratio)

temp_final$fishers_q = as.numeric(p.adjust(temp_final$fishers_exact, method = "fdr"))

# odds ratio table
temp_final_odds <- reshape2::dcast(temp_final, gene1 ~ gene2, value.var="odds_ratio_log")
temp_final_odds = as.data.frame(temp_final_odds)
rownames(temp_final_odds) <- temp_final_odds$gene1
temp_final_odds$gene1 <- NULL
temp_final_odds[is.na(temp_final_odds)] <- 1
temp_final_odds = as.matrix(temp_final_odds)

# q value table to pass into p.mat
temp_final_q <- reshape2::dcast(temp_final, gene1 ~ gene2, value.var="fishers_q")
temp_final_q = as.data.frame(temp_final_q)
rownames(temp_final_q) <- temp_final_q$gene1
temp_final_q$gene1 <- NULL
temp_final_q[is.na(temp_final_q)] <- 1
temp_final_q = as.matrix(temp_final_q)


pdf(file = "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/MetaAML_mutation_correlation_de_novo_fishers_exact.pdf", width = 7.5, height = 7.5)

corrplot(temp_final_odds, is.corr = F, type="upper", order="hclust",tl.col="black", outline = F, addgrid.col = "lightgrey",
         col = brewer.pal(n = 8, name = "RdBu"), diag=FALSE, p.mat = temp_final_q, insig = "label_sig",
         sig.level = c(.001, .01, .1), pch.cex = .9, pch.col = "black", na.label = "square", na.label.col = "white")
dev.off()

write.csv(temp_final, "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/oods_ratio_and_fishers_results.csv", row.names=FALSE)

