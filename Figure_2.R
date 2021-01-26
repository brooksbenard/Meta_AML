# ========================================================================================================================================== #
  # Figure_2.R
  # Author : Brooks Benard, bbenard@stanford.edu
  # Date: 01/25/2021
  # Description: This script will perform survival and mutation co-occurence analyses as seen in Figure 2 and related suppliments of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
  # ========================================================================================================================================== #

# ============== #
# Load libraries #
# ============== #
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('tydyr')) install.packages('tydyr'); library('tydyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('cometExactTest')) install.packages('cometExactTest'); library('cometExactTest')
if (!require('discover')) install.packages('discover'); library('discover')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('maditr')) install.packages('maditr'); library('maditr')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('epitools')) install.packages('epitools'); library('epitools')
if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('rlist')) install.packages('rlist'); library('rlist')

dir.create("~/Desktop/MetaAML_results/Figure_2")
dir.create("~/Desktop/MetaAML_results/Figure_2/Supplimental")

# load from data frame created in the Figure_1 script
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

#### Individual and Pairwise genotype associations with clinical features ####
sub = subset(final_data_matrix, mut_freq_gene > 50 & final_data_matrix$Subset == "de_novo" & final_data_matrix$Gene != "MLL")

# make sure that the FLT3 symbols are annotated well
for(i in 1:nrow(sub)){
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "ITD"){
    sub$Gene[i] <- "FLT3-ITD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "SNV"){
    sub$Gene[i] <- "FLT3-TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "Deletion"){
    sub$Gene[i] <- "FLT3-TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "INDEL"){
    sub$Gene[i] <- "FLT3-ITD"
  }
}

genes = data.frame(unique(sub$Gene))

# create directories for linear regression and binary comparisions
# regression
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Age/continuous")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_WBC/continuous")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Platelet/continuous")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_LDH/continuous")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Hemoglobin/continuous")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_PB_blast_percent/continuous")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_BM_blast_percent/continuous")

# binary
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Age/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_WBC/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Platelet/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_LDH/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Hemoglobin/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_PB_blast_percent/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_BM_blast_percent/discrete")

# summary plots for each variable ####
sub$VAF <- as.numeric(sub$VAF)
sub$Age <- as.numeric(sub$Age)
sub$Hemoglobin <- as.numeric(sub$Hemoglobin)
sub$Platelet <- as.numeric(sub$Platelet)
sub$WBC <- as.numeric(sub$WBC)
sub$LDH <- as.numeric(sub$LDH)
sub$BM_blast_percent <- as.numeric(sub$BM_blast_percent)
sub$PB_blast_percent <- as.numeric(sub$PB_blast_percent)

sub$Gene <- as.factor(sub$Gene)

variables = c("Age", "Platelet", "Hemoglobin", "WBC", "LDH", "BM_blast_percent", "PB_blast_percent")

colors = c("#bdbdbd",  "#fb6a4a", "#fe9929", "#d4b9da", "#238443", "#e6ab02", "#8dd3c7")

# age
ggplot(sub, aes(x=reorder(Gene, -Age, median, na.rm = T), y=Age)) +
  geom_boxplot(fill="#bdbdbd", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("Age") +
  geom_hline(yintercept=mean(sub$Age, na.rm = T), linetype="dashed", color = "#99000d")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Age_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 

# Platelet
ggplot(sub, aes(x=reorder(Gene, -Platelet, median, na.rm = T), y=Platelet)) +
  geom_boxplot(fill="#fb6a4a", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("Platelet (1e3/μL)") +
  geom_hline(yintercept=median(sub$Platelet, na.rm = T), linetype="dashed", color = "#99000d")+
  facet_zoom(ylim = c(0, 200))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Platelet_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 

# Hemoglobin
ggplot(sub, aes(x=reorder(Gene, -Hemoglobin, median, na.rm = T), y=Hemoglobin)) +
  geom_boxplot(fill="#fe9929", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("Hemoglobin (g/dL)") +
  geom_hline(yintercept=mean(sub$Hemoglobin, na.rm = T), linetype="dashed", color = "#99000d") +
  facet_zoom(ylim = c(0, 18))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Hemoglobin_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 

# WBC
ggplot(sub, aes(x=reorder(Gene, -WBC, median, na.rm = T), y=WBC)) +
  geom_boxplot(fill="#d4b9da", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("WBC") +
  geom_hline(yintercept=mean(sub$WBC, na.rm = T), linetype="dashed", color = "#99000d")+
  facet_zoom(ylim = c(0, 125))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/WBC_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 

# LDH
ggplot(sub, aes(x=reorder(Gene, -LDH, median, na.rm = T), y=LDH)) +
  geom_boxplot(fill="#238443", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("LDH") +
  geom_hline(yintercept=mean(sub$LDH, na.rm = T), linetype="dashed", color = "#99000d")+
  facet_zoom(ylim = c(0, 1500))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/LDH_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 

# BM_blast_percent
ggplot(sub, aes(x=reorder(Gene, -BM_blast_percent, median, na.rm = T), y=BM_blast_percent)) +
  geom_boxplot(fill="#e6ab02", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("BM Blast %") +
  geom_hline(yintercept=mean(sub$BM_blast_percent, na.rm = T), linetype="dashed", color = "#99000d") 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/BM_blast_percent_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 

# PB_blast_percent
ggplot(sub, aes(x=reorder(Gene, -PB_blast_percent, median, na.rm = T), y=PB_blast_percent)) +
  geom_boxplot(fill="#8dd3c7", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("PB Blast %") +
  geom_hline(yintercept=mean(sub$PB_blast_percent, na.rm = T), linetype="dashed", color = "#99000d")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/PB_blast_percent_vs_gene.png", dpi = 300, width = 10, height = 5, units = "in") 




# pairwise genotype correlations with clinical features ####
gene_pairs=as.data.frame(t(combn(as.vector(unique(sub$Gene)),2)))

sub2 = select(sub, Sample, Gene, PatientId, Age, Platelet, Hemoglobin, WBC, LDH, BM_blast_percent, PB_blast_percent)

gene_pairs_list = list()
z = 1

for(i in 1:length(variables)){
  print(i)
  var = as.character(variables[i])
  print(var)
  sub3 = sub2
  
  names(sub3)[names(sub3) == variables[i]] <- "Variable"
  
  for(j in 1:nrow(gene_pairs)){
    # print(j)
    gene_1 = subset(sub3, sub3$Gene == as.character(gene_pairs[j,1]))
    gene_2 = subset(sub3, sub3$Gene == as.character(gene_pairs[j,2]))
    
    overlap = inner_join(gene_1, gene_2, by = "Sample")
    overlap = overlap$Sample
    
    temp = sub3
    
    temp$genotype = ifelse(temp$Sample %in% overlap, "Double", "Other")
    
    temp = unique(select(temp, Sample, Variable, genotype))
    
    n = as.numeric(length(unique(temp$Sample)))
    
    n_check = as.numeric(length(which(temp$genotype == "Double")))
    
    if(n_check >= 5){
      print(j)
      # calculate a p-value and effect size for the difference in clinical feaures based on genotype differences
      p_val =  wilcox.test(as.numeric(temp$Variable) ~ temp$genotype, alternative = "two.sided")$p.value
      
      effect_size = cohens_d(as.numeric(Variable) ~ genotype, data = temp)
      n_n1 = as.numeric(length(which(temp$genotype == "Double")))
      n_n2 = as.numeric(length(which(temp$genotype == "Other")))
      effect_size_ci = cohen.d.ci(d = effect_size, n = n, n1 = n_n1, n2 = n_n2)
      
      
      variable_stats <- data.frame(matrix(NA, nrow = 1, ncol = 8))
      names(variable_stats) <- c("Variable", "Gene_1", "Gene_2", "effect_size", "CI_lower", "CI_upper", "p_value", "n_double")
      
      variable_stats[1,1] <- paste(var)
      variable_stats[1,2] <- as.character(gene_pairs[j,1])
      variable_stats[1,3] <- as.character(gene_pairs[j,2])
      variable_stats[1,4] <- effect_size
      variable_stats[1,5] <- effect_size_ci[1]
      variable_stats[1,6] <- effect_size_ci[3]
      variable_stats[1,7] <- p_val
      variable_stats[1,8] <- n_check
      
      # Add each list in the loop to a list of lists
      gene_pairs_list[[z]] = variable_stats 
      
      z = z + 1
    } 
  }
}

gene_pairs_list_final = do.call(rbind, gene_pairs_list)

# correct for multiple hupotheses by variable
var1_adj = list()

z = 1

for(i in 1:length(variables)){
  print(i)
  variable = as.character(variables[i])
  
  var1 = subset(gene_pairs_list_final, gene_pairs_list_final$Variable == variable)
  
  var1$q_val = p.adjust(var1$p_value)
  
  var1_adj[[z]] = var1 
  
  z = z + 1
  
}
gene_pairs_list_final_adj <- do.call(rbind, var1_adj)

# annotating color for the significant points
gene_pairs_list_final_adj$sig_color = "not_sig_ES"

for(i in 1:nrow(gene_pairs_list_final_adj)){
  if(gene_pairs_list_final_adj$q_val[i] <= 0.05 & gene_pairs_list_final_adj$effect_size[i] > 0){
    gene_pairs_list_final_adj$sig_color[i] = "sig_pos_ES"
  }
  if(gene_pairs_list_final_adj$q_val[i] <= 0.05 & gene_pairs_list_final_adj$effect_size[i] < 0){
    gene_pairs_list_final_adj$sig_color[i] = "sig_neg_ES"
  }
}


gene_pairs_list_final_adj$point_label = paste(gene_pairs_list_final_adj$Gene_1, " + ", gene_pairs_list_final_adj$Gene_2, sep = "")

for(i in 1:nrow(gene_pairs_list_final_adj)){
  if(gene_pairs_list_final_adj$q_val[i] > 0.05){
    gene_pairs_list_final_adj$point_label[i] = ""
  }
}

gene_pairs_list_final_adj$sig_color = as.factor(gene_pairs_list_final_adj$sig_color)

max_p = max(-log10(gene_pairs_list_final_adj$p_value))

gene_pairs_list_final_adj$Variable = factor(gene_pairs_list_final_adj$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM_blast_percent', 'PB_blast_percent', 'Age'))

# plot results per variable
p = ggplot(gene_pairs_list_final_adj, aes(x=gene_pairs_list_final_adj$effect_size, y=-log10(gene_pairs_list_final_adj$p_value), color = sig_color, size = n_double)) +
  # geom_label_repel(aes(label=point_label),box.padding = .5, hjust=0, vjust=1.5, size = 2.5) +
  geom_hline(yintercept = 3.563053,  linetype = "dashed", color = "lightgrey") +
  geom_point(alpha = 0.75) +
  geom_point(shape = 21, color = "black", alpha = 0.25)  +
  scale_size_area(max_size = 10,breaks=c(10,25,50,100,250,350)) +
  geom_label_repel(aes(label=point_label),box.padding = .5, hjust=1, vjust=1.5, size = 2.5) +
  # geom_hline(yintercept = 3.563053,  linetype = "dashed", color = "lightgrey") +
  theme_cowplot() +
  scale_colour_manual(values = c("sig_pos_ES"="#b35806", "sig_neg_ES" = "#542788", "not_sig_ES"="lightgrey")) +
  theme(legend.position="right") +
  # geom_hline(yintercept = 1.30103,  linetype = "dashed", color = "darkgrey") +
  ylab(label= "-log10(p-value)") +
  xlab(label= "Effect Size (co-mut vs. others)") +
  # scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(color="black", size=20)) 

g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p + facet_wrap(. ~ Variable, ncol = 7) + guides(size = g, color = FALSE) +
  theme(
    # legend.title = element_blank(), 
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid")) 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Pairwise_genotype_clinical_correlations.png", dpi = 300, width = 15, height = 5, units = "in")



#### statistical analysis of mutation co-occurence ####
# subset analysis to the de novo cohort
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
dir.create("~/Desktop/MetaAML_results/Figure_2/Supplimental")
write.csv(DISCOVER_final_2, "~/Desktop/MetaAML_results/Figure_2/Supplimental/DISCOVER_results.csv", row.names=FALSE)



temp_final_q = select(DISCOVER_final_2, gene1, gene2, q.value)
temp_final_q <- reshape2::dcast(temp_final_q, gene1 ~ gene2, value.var="q.value")
temp_final_q = as.data.frame(temp_final_q)
rownames(temp_final_q) <- temp_final_q$gene1
temp_final_q$gene1 <- NULL
temp_final_q[is.na(temp_final_q)] <- 1
temp_final_q = as.matrix(temp_final_q)

pdf(file = "~/Desktop/MetaAML_results/Figure_2/Supplimental/MetaAML_mutation_correlation_de_novo.pdf", width = 7.5, height = 7.5)

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

dir.create("~/Desktop/MetaAML_results/Figure_2/Fishers")

pdf(file = "~/Desktop/MetaAML_results/Figure_2/Fishers/MetaAML_mutation_correlation_de_novo_fishers_exact.pdf", width = 7.5, height = 7.5)

corrplot(temp_final_odds, is.corr = F, type="upper", order="hclust",tl.col="black", outline = F, addgrid.col = "lightgrey",
         col = brewer.pal(n = 8, name = "RdBu"), diag=FALSE, p.mat = temp_final_q, insig = "label_sig",
         sig.level = c(.001, .01, .1), pch.cex = .9, pch.col = "black", na.label = "square", na.label.col = "white")
dev.off()

write.csv(temp_final, "~/Desktop/MetaAML_results/Figure_2/Fishers/odds_ratio_and_fishers_results.csv", row.names=FALSE)

# volcano plot ####

df = as.data.frame(temp_final)

# remove duplicate rows
df = df %>% distinct(fishers_q, .keep_all = TRUE)

df$significant = "NS"

# annotate significance
for(i in 1:nrow(df)){
  if(df$fishers_q[i] < 0.05 & df$odds_ratio[i] < 1){
    df$significant[i] = "Mutually exclusive"
  }
  if(df$fishers_q[i] < 0.05 & df$odds_ratio[i] > 1){
    df$significant[i] = "Co-occuring"
  }
}

# set labels
df$labels = NA

for(i in 1:nrow(df)){
  if(-log(df$fishers_exact[i]) > 25){
    df$labels[i] = paste(df$gene1[i], " + ", df$gene2[i])
  }
}

ggplot(df, aes(x=log(df$odds_ratio), y=-log10(df$fishers_q), color = factor(significant), size = )) +
  geom_point(alpha = .75) +
  theme_cowplot() +
  geom_label_repel(aes(label=labels),hjust=0, vjust=1.5, size = 2) +
  scale_colour_manual(values = c("Co-occuring"= "#b2182b", "Mutually exclusive"="#2166ac", "NS"="lightgrey")) + 
  theme(legend.position="none") +
  ylab(label= "-log10(q-value)") +
  xlab(label= "log(Odds Ratio)") +
  labs(title = NULL) 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/volcano_plot_co_occurence.pdf", dpi = 300, width = 3.5, height =  5, units = "in")


# co-occurence and survival ####
# MetaAML_survival_by_co_occurence
# Brooks Benard
# bbenard@stanford.edu
# 04.03.2020
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis'); library('survivalAnalysis')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('vegan')) install.packages('vegan'); library('vegan')
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')

library(corrplot)
library(RColorBrewer)

# load data and subset to desired cohort
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_sub = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_sub = subset(final_data_matrix_sub, final_data_matrix_sub$mut_freq_gene >= 50 | final_data_matrix_sub$Gene == "ETV6")
final_data_matrix_sub$Time_to_OS = (final_data_matrix_sub$Time_to_OS/365)

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(final_data_matrix_sub)){
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] == "ITD"){
    final_data_matrix_sub$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] == "SNV"){
    final_data_matrix_sub$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] == "Deletion"){
    final_data_matrix_sub$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] == "INDEL"){
    final_data_matrix_sub$Gene[i] <- "FLT3-ITD"
  }
}

# select informative columns
final_data_matrix_sub = select(final_data_matrix_sub, Sample, Gene, Time_to_OS, Censor, Cohort, mut_freq_gene)


# pairwise mutations and survival ####
# find all unique pairwise combinations of mutations
genes=as.data.frame(t(combn(as.vector(unique(final_data_matrix_sub$Gene)),2)))
genes <- genes[order(genes$V1),]

n=1
results_list = list()

plot_list = list()

for(i in 1:nrow(genes)){
  print(i)
  gene1 <- as.character(genes[i,1])
  gene2 <- as.character(genes[i,2])
  
  # find patients with mutations in the genes
  sub1 <- subset(final_data_matrix_sub, Gene == gene1)
  sub2 <- subset(final_data_matrix_sub, Gene == gene2)
  sub3 <- setDT(sub1)[(Sample) %chin% sub2$Sample]
  
  # create strata columns
  sub1$STRATA1 <- 1
  sub2$STRATA2 <- 2
  sub3$STRATA3 <- 3
  
  # subset to informative columns
  sub1 <- unique(select(sub1, Sample, STRATA1))
  sub2 <- unique(select(sub2, Sample, STRATA2))
  sub3 <- unique(select(sub3, Sample, STRATA3))
  
  # add strata to the larger data matrix
  temp_sub <- left_join(final_data_matrix_sub, sub1, by = "Sample")
  temp_sub <- left_join(temp_sub, sub2, by = "Sample")
  temp_sub <- left_join(temp_sub, sub3, by = "Sample")
  
  # create consensus strate column
  temp_sub$STRATA <- rowSums(temp_sub[,c("STRATA1", "STRATA2", "STRATA3")], na.rm=TRUE)
  
  # filter to unique patients
  temp_sub_final <- temp_sub[!duplicated(temp_sub[1]),]
  
  temp_sub_final$Time_to_OS = as.numeric(temp_sub_final$Time_to_OS)
  temp_sub_final$Censor = as.numeric(temp_sub_final$Censor)
  
  nstrata <- as.numeric(length(unique(temp_sub_final$STRATA)))
  n_pts_both <- as.numeric(length(which(temp_sub_final$STRATA == 6)))
  
  if(nstrata > 3 & n_pts_both >= 10){
    
    # perform analysis between the different groups to get individual and pairwise hazard ratios
    # co-mutated vs others
    
    dat_6 = temp_sub_final
    
    dat_6$STRATA_all = ifelse(dat_6$STRATA == 6, 6,0)
    
    model_5 <- coxph( Surv(Time_to_OS, Censor) ~ STRATA_all,
                      data = dat_6)
    forest_5=cox_as_data_frame(coxphsummary = model_5, unmangle_dict = NULL,
                               factor_id_sep = ":", sort_by = NULL)
    forest_5$gene_1 = gene1
    forest_5$gene_2 = gene2
    forest_5$gene_tested = paste(gene1, gene2, "vs_others", sep = "_")
    
    # extract the log-rank p-value for the individual comparisons
    forest_5$log_rank_p = summary(model_5)$sctest[3]
    
    # populated the results into a dataframe and apend that dataframe 
    # forest_plot_data = rbindlist(list(forest_3, forest_4, forest_5))
    
    # forest_plot_data$iteration = i
    
    # add the number oco-occuring cases
    forest_5$num_pts = n_pts_both
    
    forest_plot_data = forest_5
    
    results_list[[n]] <- forest_plot_data
    n=n+1   
    
    # if(forest_plot_data$log_rank_p[1] <= 0.05 | forest_plot_data$log_rank_p[4] <= 0.05){
    
    #___________________#
    # plot all pairwise cases where there is a significant difference from WT
    if(summary(model_5)$sctest[3] < 0.05){
      
      all = temp_sub_final
      
      
      all$OS <- with(all, Surv(Time_to_OS, Censor == 1))
      #___________________#
      
      OS <- survfit(OS ~ STRATA, data = all, conf.type = "log-log")
      
      c <- c("lightgrey", "#e08214", "#8073ac", "#B24745FF")
      l_labs <- c("neither", gene1, gene2, "both")
      
      # plots the survival
      surv_plot = ggsurvplot(OS,
                             data = all,
                             log = (OS),
                             log.rank.weights = c("survdiff"),
                             pval = F,
                             test.for.trend = F,
                             pval.method.size = 3,
                             pval.coord = c(0, 0),
                             conf.int = F,
                             censor = T,
                             surv.median.line = "none",
                             risk.table = F,
                             risk.table.title = "",
                             risk.table.fontsize = 4,
                             # risk.table.height = ,
                             risk.table.y.text = F,
                             break.time.by = 5,
                             risk.table.pos = c("out"),
                             palette = c,
                             title = paste(gene1, " and ", gene2, sep = ""),
                             xlab = "Years",
                             ylim = c(0, 1.0),
                             ylab =  "Survival Probability",
                             font.main = c(12, "plain", "#252525"),
                             pval.size = 4,
                             font.x = c(12, "plain", "#252525"),
                             font.y =  c(12, "plain", "#252525"),
                             font.legend = c(12, "plain"),
                             font.tickslab = c(12, "plain", "#252525"),
                             legend.labs = l_labs,
                             legend.title = "Mutation status",
                             legend = "none",
                             ggtheme = theme_cowplot())
      
      plot_list[[i]] = surv_plot
      
      # print(surv_plot)
      # 
      #   png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/survival_by_co_occurence/",gene1, "_", gene2, "_co_occurence_survival.png", sep = ""), res = 300, width = 6, height = 4.5, units = "in")
      # 
      #   surv_plot
      #   print(surv_plot)
      #   dev.off()
      }
    }
  }


# bind all data
temp_final_hr = as.data.frame(do.call(rbind, results_list))

# correct for mulitple hypothesis testing
temp_final_hr$q_value <- p.adjust(temp_final_hr$log_rank_p, method = "fdr")

temp_final_hr[,1:3] = NULL
temp_final_hr = temp_final_hr %>%
  select(gene_1, gene_2, gene_tested, everything())

write.csv(temp_final_hr, "~/Desktop/MetaAML_results/Data/Tables/pairwise_mutations_hazard_ratio.csv")

# summary plot of all pairwise survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = TRUE,
                            ncol = 5, nrow = 8)
ggsave("~/Desktop/MetaAML_results/Figure_2/Supplimental/pairwise_survival_grid.pdf", plots, width = 17, height = 20)



# individual HRs ####
# list of unique individual mutations
genes_uniqe = sort(unique(final_data_matrix_sub$Gene))

# make list for populating all the survival curves
plot_list = list()

n=1
results_list = list()

for(i in genes_uniqe){
  # print(i)
  sub1 <- subset(final_data_matrix_sub, Gene == i)
  # create strata columns
  sub1$STRATA <- 1
  # subset to informative columns
  sub1 <- unique(select(sub1, Sample, STRATA))
  # add strata to the larger data matrix
  temp_sub <- left_join(final_data_matrix_sub, sub1, by = "Sample")
  # filter to unique patients
  temp_sub_final <- temp_sub[!duplicated(temp_sub[1]),]
  
  temp_sub_final$STRATA[is.na(temp_sub_final$STRATA)] <- 0
  
  temp_sub_final$Time_to_OS = as.numeric(temp_sub_final$Time_to_OS)
  temp_sub_final$Censor = as.numeric(temp_sub_final$Censor)
  
  model <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
                  data = temp_sub_final)
  forest=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                           factor_id_sep = ":", sort_by = NULL)
  forest$gene_1 = i
  forest$gene_2 = i
  forest$gene_tested = paste(i, "vs_WT", sep = "_")
  
  # extract the log-rank p-value for the individual comparisons
  forest$log_rank_p = summary(model)$sctest[3]
  
  results_list[[n]] <- forest
  n=n+1   
  
  # get the survival data ready
  temp_sub_final$OS <- with(temp_sub_final, Surv(Time_to_OS, Censor == 1))
  
  OS <- survfit(OS ~ STRATA, data = temp_sub_final, conf.type = "log-log")
  
  # extract the p-value and hazard ratio for the individual interactions
  p = round(as.numeric(summary(model)$sctest[3]), 3)
  
  p = ifelse(p < 0.001, paste0("p < 0.001"), paste("p =", p))
  
  hr = paste("HR = ", round(as.numeric(forest$HR), 2), " (", round(as.numeric(forest$Lower_CI), 2), "-", round(as.numeric(forest$Upper_CI), 2), ")", sep = "")
  
  p_hr = paste(p, "; ", hr, sep = "")
  
  c <- c("lightgrey", "#B24745FF")
  l_labs <- c("WT", i)
  
  # plots the survival
  surv_plot = ggsurvplot(OS,
                         data = temp_sub_final,
                         log = (OS),
                         log.rank.weights = c("survdiff"),
                         pval = paste0(p_hr),
                         test.for.trend = F,
                         pval.method.size = 1,
                         pval.coord = c(0,0),
                         conf.int = F,
                         censor = T,
                         surv.median.line = "none",
                         risk.table = F,
                         risk.table.title = "",
                         risk.table.fontsize = 4,
                         # risk.table.height = ,
                         risk.table.y.text = F,
                         break.time.by = 5,
                         risk.table.pos = c("out"),
                         palette = c,
                         title = paste(i),
                         xlab = "Years",
                         ylim = c(0, 1.0),
                         ylab =  "Survival Probability",
                         font.main = c(12, "plain", "#252525"),
                         pval.size = 4,
                         font.x = c(12, "plain", "#252525"),
                         font.y =  c(12, "plain", "#252525"),
                         font.legend = c(12, "plain"),
                         font.tickslab = c(12, "plain", "#252525"),
                         legend.labs = l_labs,
                         legend.title = "Mutation status",
                         legend = "none",
                         ggtheme = theme(plot.title = element_text(hjust = 0.5)))
  
  plot_list[[i]] = surv_plot
  # print(surv_plot)
  # 
  # png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_co_occurence/indiviaual_mutations/",i,"_survival.png", sep = ""), res = 300, width = 5, height = 4, units = "in")
  # 
  # surv_plot
  # print(surv_plot)
  # dev.off()
}
# bind all data
temp_final_hr_2 = as.data.frame(do.call(rbind, results_list))

# correct for mulitple hypothesis testing
temp_final_hr_2$q_value <- p.adjust(temp_final_hr_2$log_rank_p, method = "fdr")

# summary plot of all survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = TRUE,
                            ncol = 5, nrow = 6)
ggsave("~/Desktop/MetaAML_results/Figure_2/Supplimental/survival_grid.pdf", plots, width = 16, height = 17)


# plots ####
# corrplot ####
# correlation plot of individual and pairwise mutations and hazard ratios

# subset to useful data for corrplot
pairwise_data = select(temp_final_hr, HR, gene_1, gene_2, num_pts, q_value)

individual_data = select(temp_final_hr_2, HR, gene_1, gene_2, q_value)

pairwise_data = rbind.fill(pairwise_data, individual_data)

# make dataframe to populate with results
freq_genes = unique(final_data_matrix_sub$Gene)
freq_genes_matrix_hr <- data.frame(matrix(NA, nrow = length(freq_genes), ncol = length(freq_genes))) 
colnames(freq_genes_matrix_hr) = freq_genes
rownames(freq_genes_matrix_hr) = freq_genes

freq_genes_matrix_q <- data.frame(matrix(NA, nrow = length(freq_genes), ncol = length(freq_genes))) 
colnames(freq_genes_matrix_q) = freq_genes
rownames(freq_genes_matrix_q) = freq_genes

# loop through the pairwise survival results and populate the pairwise results that met criteria
for(i in 1:nrow(pairwise_data)){
  gene1 = pairwise_data$gene_1[i]
  gene2 = pairwise_data$gene_2[i]
  
  # populate HR
  g1 = which(colnames(freq_genes_matrix_hr) == gene1)
  g2 = which(colnames(freq_genes_matrix_hr) == gene2)
  
  freq_genes_matrix_hr[g1,g2] = -log(pairwise_data$HR[i])
  freq_genes_matrix_hr[g2,g1] = -log(pairwise_data$HR[i])
  
  # populate q val
  g1 = which(colnames(freq_genes_matrix_q) == gene1)
  g2 = which(colnames(freq_genes_matrix_q) == gene2)
  
  freq_genes_matrix_q[g1,g2] = pairwise_data$q_value[i]
  freq_genes_matrix_q[g2,g1] = pairwise_data$q_value[i]
  
}


freq_genes_matrix_hr[is.na(freq_genes_matrix_hr)] <- 0
freq_genes_matrix_hr = as.matrix(freq_genes_matrix_hr)

# freq_genes_matrix_q[is.na(freq_genes_matrix_q)] <- 0
freq_genes_matrix_q = as.matrix(freq_genes_matrix_q)

pdf(file = "~/Desktop/MetaAML_results/Figure_2/Supplimental/pairwise_survival_correlation_de_novo.pdf", width = 7.5, height = 7.5)

corrplot(freq_genes_matrix_hr,
         is.corr = F, 
         type="upper", 
         order="alphabet",
         tl.col="black", 
         outline = F, 
         addgrid.col = "lightgrey",
         col = brewer.pal(n = 8, name = "PRGn"),
         diag=F, 
         p.mat = freq_genes_matrix_q, 
         insig = "label_sig",
         sig.level = c(.001, .01, 0.1),
         pch.cex = .9,
         pch.col = "black", 
         na.label = "square",
         na.label.col = "white")
dev.off()



# forrest plot ####
# individual mutation's HRs

temp_final_hr_2$genes <- factor(temp_final_hr_2$gene_1, levels = temp_final_hr_2$gene_1[order(temp_final_hr_2$HR)])
temp_final_hr_2$log_rank_p = round(temp_final_hr_2$log_rank_p, 3)
temp_final_hr_2$p_text = NA

for(i in 1:nrow(temp_final_hr_2)){
  if(temp_final_hr_2$log_rank_p[i] <= 0.05 & temp_final_hr_2$log_rank_p[i] >= 0.001){
    temp_final_hr_2$p_text[i] = paste("p = ", temp_final_hr_2$log_rank_p[i], sep = "")
  } 
  if(temp_final_hr_2$log_rank_p[i] == 0){
    temp_final_hr_2$p_text[i] = paste("p < 0.001")
  } 
  if(temp_final_hr_2$log_rank_p[i] > 0.05){
    temp_final_hr_2$p_text[i] = ""
  }
}

# add functional category to the mutations for visualization purposes
temp_final_hr_2$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD", "JAK2")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33", "MLL")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(temp_final_hr_2)){
  if(temp_final_hr_2$genes[i] %in% DNA_methylation){
    temp_final_hr_2$mutation_category[i] <- "DNA Methylation"
  }
  if(temp_final_hr_2$genes[i] %in% Chromatin_cohesin){
    temp_final_hr_2$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(temp_final_hr_2$genes[i] %in% RTK_RAS_Signaling){
    temp_final_hr_2$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(temp_final_hr_2$genes[i] %in% Splicing){
    temp_final_hr_2$mutation_category[i] <- "Splicing"
  }
  if(temp_final_hr_2$genes[i] %in% Transcription){
    temp_final_hr_2$mutation_category[i] <- "Transcription"
  }
  if(temp_final_hr_2$genes[i] == "NPM1"){
    temp_final_hr_2$mutation_category[i] <- "NPM1"
  }
  if(temp_final_hr_2$genes[i] %in% Tumor_suppressors){
    temp_final_hr_2$mutation_category[i] <- "Tumor suppressors"
  }
}

# color coded by mutation category plot
ggplot(temp_final_hr_2, aes(x = reorder(genes, -HR), y = HR, label = temp_final_hr_2$p_text)) +
  geom_hline(yintercept=.5, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=3, linetype="dashed", color = "#d9d9d9") +
  geom_text(aes(genes, Upper_CI), hjust = 0, nudge_y = 0.35, size = 3) +
  ylim(0,6) +
  geom_pointrange(size = .75, stat = "identity", shape = 15, 
                  aes(x = genes, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = mutation_category)) +
  scale_color_manual(name = "", values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF")) +
  ylab("Hazard Ratio")+
  theme(legend.position = "right",
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 

# red vs grey color for significance
# temp_final_hr_2$color = ifelse(temp_final_hr_2$q_value <= 0.1, "q < 0.1", "q > 0.1")
temp_final_hr_2$color = "NS"
for(i in 1:nrow(temp_final_hr_2)){
  if(temp_final_hr_2$HR[i] > 1 & temp_final_hr_2$q_value[i] <= 0.1){
    temp_final_hr_2$color[i] = "Better"
  }
  if(temp_final_hr_2$HR[i] < 1 & temp_final_hr_2$q_value[i] <= 0.1){
    temp_final_hr_2$color[i] = "Worse"
  }
}

ggplot(temp_final_hr_2, aes(x = reorder(genes, -HR), y = HR, label = temp_final_hr_2$p_text)) +
  geom_hline(yintercept=.5, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=3, linetype="dashed", color = "#d9d9d9") +
  geom_text(aes(genes, Upper_CI), hjust = 0, nudge_y = 0.35, size = 3) +
  ylim(0,5.5) +
  geom_pointrange(size = .75, stat = "identity", shape = 15,
                  aes(x = genes, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = as.factor(color))) +
  scale_color_manual(name = "", values = c("Worse" = "#1b7837", "Better" = "#762a83", "NS" = "#737373")) +
  ylab("Hazard Ratio")+
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/forrest_plot_de_novo.pdf", dpi = 300, width = 7.5, height = 7.5, units = "in")




# volcano plot ####

df = as.data.frame(temp_final_hr)

# remove duplicate rows
df = df %>% distinct(df$q_value, .keep_all = TRUE)

df$significant = "NS"

# annotate significance
for(i in 1:nrow(df)){
  if(df$q_value[i] < 0.05 & df$HR[i] < 1){
    df$significant[i] = "Favorable"
  }
  if(df$q_value[i] < 0.05 & df$HR[i] > 1){
    df$significant[i] = "Unfavorable"
  }
}

# set labels
df$labels = NA

for(i in 1:nrow(df)){
  if(df$q_value[i] < 0.05){
    df$labels[i] = paste(df$gene_1[i], " + ", df$gene_2[i])
  }
}

# tp53 is a huge outlier so for the plot remove it and then later add it ot the upper right corner
# df = df[!(df$gene_1=="TP53" & df$gene_2=="TP53"),]

p = ggplot(df, aes(x=df$HR, y=-log10(df$q_value), color = factor(significant), size = num_pts)) +
  geom_point(alpha = .75) +
  geom_label_repel(aes(label=labels),hjust=0, vjust=1.5, size = 2) +
  scale_colour_manual(values = c("Favorable"= "#1b7837", "Unfavorable" = "#762a83", "NS"="#737373")) + 
  theme_cowplot() +
   theme(legend.position="nont") +
  ylab(label= "-log10(q-value)") +
  xlab(label= "Hazard Ratio") +
  labs(title = NULL) +
  theme(plot.title = element_text(color="black", size=20))

g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p  + guides(size = g, color = FALSE) 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/volcano_plot_co_occurence_survival.pdf", dpi = 300, width = 3.5, height =  5, units = "in")



# scattterplot of HR and odds ratio ####
# read in the results from MetaAML_co_occurence_mutual_exlclusive_analysis.R
odds_ratio = read.csv("~/Desktop/MetaAML_results/Figure_2/oods_ratio_and_fishers_results.csv")
odds_ratio = select(odds_ratio, gene1, gene2, odds_ratio, fishers_q)
colnames(odds_ratio) = c("gene_1", "gene_2", "odds_ratio", "q_value_odds")

hazard_ratio = pairwise_data[,c(1:3, 5)]
colnames(hazard_ratio)[4] = "q_value_HR"

hr_odds = inner_join(odds_ratio, hazard_ratio, by = c("gene_1", "gene_2"))

# hr_odds$color = as.factor(ifelse(hr_odds$q_value_HR < 0.1, "HR q < 0.1","HR q > 0.1"))

hr_odds$color = "NS"

for(i in 1:nrow(hr_odds)){
  if(hr_odds$q_value_HR[i] < 0.05 & hr_odds$HR[i] > 1){
    hr_odds$color[i] = "HR > 1\nq < 0.05\n"
  }
  if(hr_odds$q_value_HR[i] < 0.05 & hr_odds$HR[i] < 1){
    hr_odds$color[i] = "HR < 1\nq < 0.05\n"
  }
}

# add annotations for most significant interactions
hr_odds$labels = NA

for(i in 1:nrow(hr_odds)){
  if(hr_odds$q_value_HR[i] < 0.05){
    if(log(hr_odds$HR[i]) > 0.18 | log(hr_odds$HR[i]) < -0.2 | log(hr_odds$odds_ratio[i]) < -1 | log(hr_odds$odds_ratio[i]) > 1){
      hr_odds$labels[i] = paste(hr_odds$gene_1[i], "+",  hr_odds$gene_2[i])
    } 
  }
}

ggplot(hr_odds, aes(x = log(odds_ratio), y = log(HR), color = as.factor(color))) +
  geom_point(aes(color = color, shape = ), size = 3, alpha = 0.75) +
  xlab("log(Odds Ratio)") +
  ylab("log(Hazard Ratio)") +
  scale_color_manual(values = c("#1b7837",  "#762a83", "grey")) +
  geom_point(shape = 1, size =  3,colour = "black") +
  geom_hline(yintercept=log(1), linetype="dashed", color = "grey") +
  geom_vline(xintercept=log(1), linetype="dashed", color = "grey") +
  geom_label_repel(aes(label=labels),
                   size = 2.5) +
  theme_cowplot() +
  theme(
    legend.position = "none")
  


ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/correlation_of_odds_ratio_and_HR_de_novo.pdf", dpi = 300, width = 5, height = 5, units = "in")

write.csv(hr_odds, "~/Desktop/MetaAML_results/Data/Tables/correlation_of_odds_ratio_and_HR_de_novo.pdf.csv")

# test the enrichment of poor outcomes based on odds ratio

# for odds ratio greater than 1
one = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) > 0 & hr_odds$q_value_odds < .1)))
two = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) > 0)))
three = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & hr_odds$q_value_odds < .1)))
four = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) < 0)))

or_hr <- data.frame(matrix(NA, nrow = 2, ncol = 2))
names(or_hr) <- c("or", "hr")

or_hr[1,1] <- 12
or_hr[1,2] <- 9
or_hr[2,1] <- 44
or_hr[2,2] <- 29

prop.test(table(or_hr$or, or_hr$hr), correct=FALSE) 


# for odds ratio less than 1
one = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) > 0 & hr_odds$q_value_odds < .1)))
two = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) > 0)))
three = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & hr_odds$q_value_odds < .1)))
four = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) < 0)))

or_hr <- data.frame(matrix(NA, nrow = 2, ncol = 2))
names(or_hr) <- c("or", "hr")

or_hr[1,1] <- 11
or_hr[1,2] <- 4
or_hr[2,1] <- 29
or_hr[2,2] <- 31

prop.test(table(or_hr$or, or_hr$hr), correct=FALSE) 