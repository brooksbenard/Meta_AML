# ========================================================================================================================================== #
# Figure_2.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 03/16/2021
# Description: This script will perform survival and mutation co-occurence analyses as seen in Figure 2 and related suppliments of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# ============== #
# Load libraries #
# ============== #
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('tydyr')) install.packages('tydyr', dependencies = TRUE); library('tydyr')
if (!require('dplyr')) install.packages('dplyr', dependencies = TRUE); library('dplyr')
if (!require('cometExactTest')) install.packages('cometExactTest', dependencies = TRUE); library('cometExactTest')
if (!require('discover')) install.packages('discover'); library('discover')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('maditr')) install.packages('maditr'); library('maditr')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('epitools')) install.packages('epitools'); library('epitools')
if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis', dependencies = TRUE); library('survivalAnalysis')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('vegan')) install.packages('vegan'); library('vegan')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('ggforce')) install.packages('ggforce'); library('ggforce')
if (!require('rstatix')) install.packages('rstatix'); library('rstatix')
if (!require('effsize')) install.packages('effsize'); library('effsize')
if (!require('psych')) install.packages('psych'); library('psych')

dir.create("~/Desktop/MetaAML_results/Figure_2")
dir.create("~/Desktop/MetaAML_results/Figure_2/Supplimental")

# load from data frame created in the Figure_1 script
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

#### Individual and Pairwise genotype associations with clinical features ####
sub = subset(final_data_matrix, mut_freq_gene > 50 & final_data_matrix$Subset == "de_novo")

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
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "other"){
    sub$Gene[i] <- "FLT3-TKD"
  }
}

genes = data.frame(unique(sub$Gene))

# create directories for linear regression and binary comparisions
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Age")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_WBC")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Platelet")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_LDH")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Hemoglobin")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_PB_blast_percent")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_BM_blast_percent")

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

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Age_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

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

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Platelet_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

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

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Hemoglobin_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

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

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/WBC_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

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

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/LDH_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

# BM_blast_percent
ggplot(sub, aes(x=reorder(Gene, -BM_blast_percent, median, na.rm = T), y=BM_blast_percent)) +
  geom_boxplot(fill="#e6ab02", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("BM Blast %") +
  geom_hline(yintercept=mean(sub$BM_blast_percent, na.rm = T), linetype="dashed", color = "#99000d") 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/BM_blast_percent_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

# PB_blast_percent
ggplot(sub, aes(x=reorder(Gene, -PB_blast_percent, median, na.rm = T), y=PB_blast_percent)) +
  geom_boxplot(fill="#8dd3c7", color="black") +
  theme_cowplot(font_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                   size = 10, hjust = .5)) +
  xlab(label = NULL) +
  ylab("PB Blast %") +
  geom_hline(yintercept=mean(sub$PB_blast_percent, na.rm = T), linetype="dashed", color = "#99000d")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/PB_blast_percent_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 




# pairwise genotype correlations with clinical features ####
gene_pairs=as.data.frame(t(combn(as.vector(unique(sub$Gene)),2)))

sub2 = dplyr::select(sub, Sample, Gene, PatientId, Age, Platelet, Hemoglobin, WBC, LDH, BM_blast_percent, PB_blast_percent)

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
    
    temp = unique(dplyr::select(temp, Sample, Variable, genotype))
    
    n = as.numeric(length(unique(temp$Sample)))
    
    n_check = as.numeric(length(which(temp$genotype == "Double")))
    
    if(n_check >= 5){
      print(j)
      # calculate a p-value and effect size for the difference in clinical feaures based on genotype differences
      p_val =  wilcox.test(as.numeric(temp$Variable) ~ temp$genotype, alternative = "two.sided")$p.value
      
      effect_size = cohens_d(Variable ~ genotype, data = temp)$effsize
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

gene_pairs_list_final_adj$Variable = factor(gene_pairs_list_final_adj$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM_blast_percent', 'PB_blast_percent', 'Age'))

# plot results per variable
p = ggplot(gene_pairs_list_final_adj, aes(x=effect_size, y=-log10(p_value), color = sig_color, size = n_double)) +
  theme_cowplot() +
  geom_hline(yintercept = 3.563053,  linetype = "dashed", color = "lightgrey") +
  geom_point(alpha = 0.75) +
  xlim(-1,2) +
  geom_point(shape = 21, color = "black", alpha = 0.25)  +
  geom_label_repel(aes(label=point_label), size = 2.5, force = 75, max.overlaps = 5) +
  scale_size_area(max_size = 10,breaks=c(10,25,50,100,250,350)) +
  scale_colour_manual(values = c("sig_pos_ES"="#b35806", "sig_neg_ES" = "#542788", "not_sig_ES"="lightgrey")) +
  theme(legend.position="right") +
  ylab(label= "-log10(p-value)") +
  xlab(label= "Effect Size (co-mut vs. others)") +
  theme(plot.title = element_text(color="black", size=20)) 

g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p = p + facet_wrap(. ~ Variable, ncol = 7) + guides(size = g, color = FALSE) +
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid")) 

ggsave(plot = p, filename = "~/Desktop/MetaAML_results/Figure_2/Pairwise_genotype_clinical_correlations.png", dpi = 300, width = 15, height = 5, units = "in")



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
# events <- discover.matrix(final_data_matrix_2_sub)
# subset <- rowSums(final_data_matrix_2_sub) > 25
# result.mutex <- pairwise.discover.test(events[subset, ])
# result.mutex
# print(result.mutex, fdr.threshold=0.05)
# result.mutex = as.data.frame(result.mutex)

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
    odds_ratio = 0.06982398

  }
  
  # store the results in a dataframe
  odds <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  names(odds) <- c("gene1", "gene2", "odds_ratio", "fishers_exact", "n_cooccur")
  
  odds[1,1] <- gene_1
  odds[1,2] <- gene_2
  odds[1,3] <- odds_ratio
  odds[1,4] <- fishers_p
  odds[1,5] <- temp$Freq[4]
  
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
    odds_ratio = 0.06982398

  }
  
  
  # store the results in a dataframe
  odds <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  names(odds) <- c("gene1", "gene2", "odds_ratio", "fishers_exact", "n_cooccur")
  
  odds[1,1] <- gene_1
  odds[1,2] <- gene_2
  odds[1,3] <- odds_ratio
  odds[1,4] <- fishers_p
  odds[1,5] <- temp$Freq[4]
  
  results_list[[n]] <- odds
  n=n+1   
}
temp_final_2 = as.data.frame(do.call(rbind, results_list))

temp_final = unique(rbind(temp_final_1, temp_final_2))
# use temp_final later in visualizing the fishers exact results

# Fisher's exact test  ####
# plot results from Fisher's Exact analysis
# log transform the odds ratio so that it can be easily visualized in the corrplot
temp_final$odds_ratio_log = -log(temp_final$odds_ratio)

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
  if(-log10(df$fishers_exact[i]) > 10){
    df$labels[i] = paste0(df$gene1[i], " + ", df$gene2[i])
  }
}

p = ggplot(df, aes(x=log(odds_ratio), y=-log10(fishers_q), color = significant, size = n_cooccur)) +
  theme_cowplot() +
  geom_point(alpha = .75) +
  geom_label_repel(aes(label=labels),nudge_y = 5, size = 3, force = 25, max.overlaps = 5) +
  scale_colour_manual(values = c("Co-occuring"= "#b2182b", "Mutually exclusive"="#2166ac", "NS"="lightgrey")) + 
  theme(legend.position = c(0.05,.85)) +
  ylab(label= "-log10(q-value)") +
  xlab(label= "log(Odds Ratio)") +
  scale_size_area(max_size = 5,breaks=c(25,50,100,200,300)) +
  labs(title = NULL) 

g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p + guides(size = g, color = FALSE) 

 ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/volcano_plot_co_occurence.pdf", dpi = 300, width = 3.5, height =  5, units = "in")


# co-occurence and survival ####

# load data and subset to desired cohort
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_sub = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_sub = subset(final_data_matrix_sub, final_data_matrix_sub$mut_freq_gene >= 100)
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
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] == "other"){
    final_data_matrix_sub$Gene[i] <- "FLT3-TKD"
  }
}


# select informative columns
final_data_matrix_sub = select(final_data_matrix_sub, Sample, Gene, Time_to_OS, Censor, Cohort)


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
  sub3 <- setDT(sub1)[(Sample) %in% sub2$Sample]
  
  n_pts = as.numeric(n_distinct(sub3$Sample))
  
  if(n_pts >= 10){
    
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
    temp_sub_final <- temp_sub[!duplicated(temp_sub$Sample),]
    
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
      
      # create new data frame to populate results
      forest_5 <- data.frame(matrix(NA, nrow = 1, ncol = 6))
      names(forest_5) <- c("gene_1", "gene_2", "gene_tested", "HR", "lower_95", "upper_95")
      
      forest_5$gene_1 = gene1
      forest_5$gene_2 = gene2
      forest_5$gene_tested = paste(gene1, gene2, "vs_others", sep = "_")
      forest_5$HR = round(exp(coef(model_5)), 2)
      forest_5$lower_95 = round(exp(confint(model_5)[1]), 2)
      forest_5$upper_95 = round(exp(confint(model_5)[2]), 2)
      
      # extract the log-rank p-value for the individual comparisons
      forest_5$log_rank_p = summary(model_5)$sctest[3]
      
      # add the number of co-occuring cases
      forest_5$num_pts = n_pts_both
      
      forest_plot_data = forest_5
      
      results_list[[n]] <- forest_plot_data
      n=n+1   
      
      # plot all pairwise cases where there is a significant difference from WT
      if(summary(model_5)$sctest[3] < 0.021){
        
        # extract the p-value and hazard ratio for the individual interactions
        p = round(as.numeric(summary(model_5)$sctest[3]), 3)
        
        p = ifelse(p < 0.001, paste0("p < 0.001"), paste("p =", p))
        
        hr = paste("HR = ", round(as.numeric(forest_5$HR), 2), " (", round(as.numeric(forest_5$lower_95), 2), "-", round(as.numeric(forest_5$upper_95), 2), ")", sep = "")
        
        p_hr = paste("WT vs. co-mut\n", p, "\n", hr, sep = "")
        
        all = dat_6
        
        OS_fit <- survfit(Surv(Time_to_OS, Censor) ~ 1, data=all)
        OS_trt_fit <- survfit(Surv(Time_to_OS, Censor) ~ STRATA, data=all, conf.type = "log-log")
        
        c <- c("lightgrey", "#e08214", "#8073ac", "#B24745FF")
        l_labs <- c("neither", gene1, gene2, "both")
        
        # plots the survival
        surv_plot = ggsurvplot(OS_trt_fit,
                               data = all,
                               log = (OS_trt_fit),
                               log.rank.weights = c("survdiff"),
                               pval = paste0(p_hr),
                               test.for.trend = F,
                               pval.method.size = .5,
                               pval.coord = c(7.5,.85),
                               conf.int = F,
                               censor = T,
                               surv.median.line = "none",
                               risk.table = F,
                               risk.table.title = "",
                               risk.table.fontsize = 4,
                               risk.table.y.text = F,
                               break.time.by = 5,
                               risk.table.pos = c("out"),
                               palette = c,
                               title = paste(gene1, " + ", gene2, sep = ""),
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
}


# bind all data
temp_final_hr = as.data.frame(do.call(rbind, results_list))

# correct for mulitple hypothesis testing
temp_final_hr$q_value <- p.adjust(temp_final_hr$log_rank_p, method = "fdr")

# temp_final_hr[,1:3] = NULL
temp_final_hr = temp_final_hr %>%
  select(gene_1, gene_2, gene_tested, everything())

write.csv(temp_final_hr, "~/Desktop/MetaAML_results/Data/Tables/pairwise_mutations_hazard_ratio.csv")

# summary plot of all pairwise survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = TRUE,
                            ncol = 4, nrow = 6)
ggsave("~/Desktop/MetaAML_results/Figure_2/pairwise_survival_grid.pdf", plots, width = 12.75, height = 15)



# individual HRs ####
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_sub = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_sub = subset(final_data_matrix_sub, final_data_matrix_sub$mut_freq_gene >= 50)
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
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] == "other"){
    final_data_matrix_sub$Gene[i] <- "FLT3-TKD"
  }
}

# list of unique individual mutations
genes_uniqe =  c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1","NRAS", "KRAS", "FLT3-ITD", "FLT3-TKD", "PTPN11", "STAG2")

final_data_matrix_sub_simple = na.omit(distinct(final_data_matrix_sub, Sample, Time_to_OS, Censor))

# make list for populating all the survival curves
plot_list = list()

n=1
results_list = list()

for(i in genes_uniqe){
  print(i)
  # find mutated pts for the given gene
  # gene = as.character(genes_uniqe[i])
  sub1 <- subset(final_data_matrix_sub, Gene == i)
  mut_pts = unique(sub1$Sample)
  
  temp = final_data_matrix_sub_simple
  temp$STRATA = 0
  
  for(j in 1:nrow(temp)){
    if(temp$Sample[j] %in% mut_pts) {temp$STRATA[j] = 1} 
  }

  temp$Time_to_OS = as.numeric(temp$Time_to_OS)
  temp$Censor = as.numeric(temp$Censor)
  
  model <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
                  data = temp)
  
  forest <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  names(forest) <- c("gene", "test_groups", "HR", "lower_95", "upper_95")
  
  forest$gene = paste(i)
  forest$test_groups = paste(i, "vs_WT", sep = "_")
  forest$HR = round(exp(coef(model)), 2)
  forest$lower_95 = round(exp(confint(model)[1]), 2)
  forest$upper_95 = round(exp(confint(model)[2]), 2)

  # extract the log-rank p-value for the individual comparisons
  forest$log_rank_p = summary(model)$sctest[3]
  
  results_list[[n]] <- forest
  n=n+1   
  
  # get the survival data ready
  OS_fit <- survfit(Surv(Time_to_OS, Censor) ~ 1, data=temp)
  OS_trt_fit <- survfit(Surv(Time_to_OS, Censor) ~ STRATA, data=temp, conf.type = "log-log")
  
  # extract the p-value and hazard ratio for the individual interactions
  p = round(as.numeric(summary(model)$sctest[3]), 3)
  
  p = ifelse(p < 0.001, paste0("p < 0.001"), paste("p =", p))
  
  hr = paste("HR = ", round(as.numeric(forest$HR), 2), " (", round(as.numeric(forest$lower_95), 2), "-", round(as.numeric(forest$upper_95), 2), ")", sep = "")
  
  p_hr = paste(p, "; ", hr, sep = "")
  
  c <- c("lightgrey", "#B24745FF")
  # l_labs <- c("WT", i)
  
  # plots the survival
  surv_plot = ggsurvplot(OS_trt_fit,
                         data = temp,
                         log = (OS_trt_fit),
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
                         legend.title = "Mutation status",
                         legend = "none",
                         ggtheme = theme_cowplot())
  
  plot_list[[i]] = surv_plot
}
# bind all data
temp_final_hr_2 = as.data.frame(do.call(rbind, results_list))

# correct for mulitple hypothesis testing
temp_final_hr_2$q_value <- p.adjust(temp_final_hr_2$log_rank_p, method = "fdr")

write.csv(temp_final_hr_2, "~/Desktop/MetaAML_results/Figure_2/Supplimental/gene_hazard_ratio.csv")

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
# pairwise_data$q_value = round(pairwise_data$q_value, 2)
# individual_data = select(temp_final_hr_2, HR, gene_1, gene_2, q_value)

# pairwise_data = rbind.fill(pairwise_data, individual_data)

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

# volcano plot ####
df = as.data.frame(temp_final_hr)

# remove duplicate rows
# df = df %>% distinct(df$q_value, .keep_all = TRUE)

df$significant = "NS"

# annotate significance
for(i in 1:nrow(df)){
  if(df$q_value[i] <= 0.1 & df$HR[i] < 1){
    df$significant[i] = "Favorable"
  }
  if(df$q_value[i] <= 0.1 & df$HR[i] > 1){
    df$significant[i] = "Unfavorable"
  }
}

# set labels
df$labels = NA

for(i in 1:nrow(df)){
  if(df$q_value[i] < 0.01 & df$HR[i] > 1.15){
    df$labels[i] = paste(df$gene_1[i], " + ", df$gene_2[i])
  }
  if(df$q_value[i] < 0.05 & df$HR[i] < 0.95){
    df$labels[i] = paste(df$gene_1[i], " + ", df$gene_2[i])
  }
}

p = ggplot(df, aes(x=df$HR, y=-log10(df$q_value), color = factor(significant), size = num_pts)) +
  theme_cowplot() +
  geom_hline(yintercept=1, linetype="dashed", color = "#d9d9d9") +
  geom_point(alpha = .75) +
  geom_label_repel(aes(label=labels),size = 3, force = 50) +
  scale_colour_manual(values = c("Favorable"= "#1b7837", "Unfavorable" = "#762a83", "NS"="#737373")) + 
  theme_cowplot() +
  theme(legend.position = c(0.05,.85)) +
  ylab(label= "-log10(q-value)") +
  xlab(label= "Hazard Ratio") +
  labs(title = NULL) +
  theme(plot.title = element_text(color="black", size=20)) 
g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p  + guides(size = g, color = FALSE) 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/volcano_plot_co_occurence_survival.pdf", dpi = 300, width = 3.5, height =  5, units = "in")



# scattterplot of HR and odds ratio ####
# read in the results from MetaAML_co_occurence_mutual_exlclusive_analysis.R
odds_ratio = read.csv("~/Desktop/MetaAML_results/Figure_2/Fishers/odds_ratio_and_fishers_results.csv")
odds_ratio = select(odds_ratio, gene1, gene2, odds_ratio, fishers_q)
colnames(odds_ratio) = c("gene_1", "gene_2", "odds_ratio", "q_value_odds")

hazard_ratio = read.csv("~/Desktop/MetaAML_results/Data/Tables/pairwise_mutations_hazard_ratio.csv")
hazard_ratio = hazard_ratio %>%
  select(gene_1, gene_2, HR, q_value)
colnames(hazard_ratio)[4] = "q_value_HR"

hr_odds = inner_join(odds_ratio, hazard_ratio, by = c("gene_1", "gene_2"))

hr_odds$color = "NS"

for(i in 1:nrow(hr_odds)){
  if(hr_odds$q_value_HR[i] < 0.1 & hr_odds$HR[i] > 1){
    hr_odds$color[i] = "HR > 1\nq < 0.1\n"
  }
  if(hr_odds$q_value_HR[i] < 0.1 & hr_odds$HR[i] < 1){
    hr_odds$color[i] = "HR < 1\nq < 0.1\n"
  }
}

# add annotations for most significant interactions
hr_odds$labels = NA

for(i in 1:nrow(hr_odds)){
  if(hr_odds$q_value_HR[i] < 0.1){
    if(log(hr_odds$HR[i]) > 0.18 | log(hr_odds$HR[i]) < -0.2 | log(hr_odds$odds_ratio[i]) < -1 | log(hr_odds$odds_ratio[i]) > 1){
      hr_odds$labels[i] = paste(hr_odds$gene_1[i], "+",  hr_odds$gene_2[i])
    } 
  }
}

ggplot(hr_odds, aes(x = log(odds_ratio), y = log(HR), color = as.factor(color))) +
  # geom_rect(data= hr_odds, inherit.aes = FALSE,
  #             aes(xmin=-Inf, xmax=+Inf, ymin=-Inf, ymax=0), 
  #             fill='#f4a582', alpha=0.5) +
  # geom_rect(data= hr_odds, inherit.aes = FALSE,
  #           aes(xmin=-Inf, xmax=+Inf, ymin=0, ymax=+Inf), 
  #           fill='#d1e5f0') +
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


ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/correlation_of_odds_ratio_and_HR_de_novo.pdf", dpi = 300, width = 4.5, height = 4.5, units = "in")

write.csv(hr_odds, "~/Desktop/MetaAML_results/Data/Tables/correlation_of_odds_ratio_and_HR_de_novo.pdf.csv")

# test the enrichment of poor outcomes based on odds ratio

# for odds ratio greater than 1
one = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) > 0 & hr_odds$q_value_odds < .1)))
two = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) > 0)))
three = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & hr_odds$q_value_odds < .1)))
four = as.numeric(nrow(subset(hr_odds, log(hr_odds$odds_ratio) > 0 & log(hr_odds$HR) < 0)))

or_hr <- data.frame(matrix(NA, nrow = 2, ncol = 2))
names(or_hr) <- c("or", "hr")

or_hr[1,1] <- one
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




# VAF Scatterplot ####
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

vaf_scatterplot_function <- function(pt_subset, gene_1_2, save_plot){
  
  # subset to desired cohorts
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix
    label1 <- "All"
  }
  if(pt_subset == "De novo"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
    label1 <- "De novo"
  }
  if(pt_subset == "Secondary"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "transformed")
    label1 <- "Secondary"
  }
  if(pt_subset == "Relapse"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
    label1 <- "Relapse"
  }
  if(pt_subset == "Therapy related"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
    label1 <- "Therapy related"
  }
  if(pt_subset == "Other"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
    label1 <- "Other"
  }
  
  # select mutations of interest
  gene_x <- gene_1_2[1]
  gene_y <- gene_1_2[2]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_x)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_y)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  gene_1_and_2=as.data.frame(gene_1_and_2) %>% distinct(Sample, VAF.x, VAF.y, .keep_all = F)
  
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF.x),]
  gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  # add point color columns for visualizing clonal/subclonal trends
  gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF.x - gene_1_and_2$VAF.y)
  
  gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
  
  # define point color
  gene_1_and_2$Clonality <- NA
  
  for(i in 1:nrow(gene_1_and_2)){
    if(gene_1_and_2$vaf_ratio[i] <= 5 & gene_1_and_2$vaf_ratio[i] >= -5){
      gene_1_and_2$Clonality[i] <- "Ambiguous"
    }
    if(gene_1_and_2$vaf_ratio[i] > 5){
      gene_1_and_2$Clonality[i] <- paste(gene_x, "first", sep = " ")
    }
    if(gene_1_and_2$vaf_ratio[i] < -5){
      gene_1_and_2$Clonality[i] <- paste(gene_y, "first", sep = " ")
    }
  }
  
  if(nrow(gene_1_and_2) > 4){
    # make the scatterplot
    
    scatter_plot = ggplot(gene_1_and_2, aes(x = gene_1_and_2$VAF.y, y = gene_1_and_2$VAF.x)) +
theme_cowplot() +
    geom_point(aes(color = Clonality), size = 3.5, alpha = 0.75) +
      xlim(0,(max(gene_1_and_2$VAF.y) + .05))+
      ylim(0,(max(gene_1_and_2$VAF.x) + .05))+
      xlab(paste(gene_y, " VAF", sep ="")) +
      ylab(paste(gene_x, " VAF", sep ="")) +  
      
      scale_color_manual(values = c("#374E55FF","#8c510a","#01665e")) +
      geom_abline(intercept = 5, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_abline(intercept = -5, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_point(shape = 1, size =  3.5,colour = "black") +
      theme(plot.title = element_text(hjust = 0.5, paste(gene_y, "vs.", gene_x, sep = "")), 
            legend.position = "right") 
    
    print(scatter_plot)
    
    if(save_plot == T){
      ggsave(filename =paste("~/Desktop/MetaAML_results/Figure_2/",gene_x, "_",gene_y,"_scatterplot.png", sep = ""), dpi = 300, width = 4.5, height = 3, units = "in")
    }
  }
}

vaf_scatterplot_function(pt_subset = "De novo", gene_1_2 = c("NRAS", "KRAS"), save_plot = T)
vaf_scatterplot_function(pt_subset = "De novo", gene_1_2 = c("NRAS", "PTPN11"), save_plot = T)



# Supplimental ####

load("~/Desktop/MetaAML_results/final_data_matrix.RData")

sub = subset(final_data_matrix, mut_freq_gene > 50 & final_data_matrix$Subset == "de_novo")

# make sure that the FLT3 symbols are the same
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
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "other"){
    sub$Gene[i] <- "FLT3-TKD"
  }
}

sub = unique(select(sub, Sample, Gene, PatientId, Age, BM_blast_percent,  PB_blast_percent, WBC, Hemoglobin, LDH, Platelet, PB_wbc_percent))

genes = data.frame(unique(sub$Gene))

variables = c("Age", "Platelet", "Hemoglobin", "WBC", "LDH", "BM_blast_percent", "PB_blast_percent")

# pairwise genotype correlations with clinical features ####
sub2 = select(sub, Sample, Gene, PatientId, Age, Platelet, Hemoglobin, WBC, LDH, BM_blast_percent, PB_blast_percent)

gene_pairs_list = list()
z = 1

for(i in 1:length(variables)){
  print(i)
  var = as.character(variables[i])
  print(var)
  sub3 = sub2
  
  names(sub3)[names(sub3) == variables[i]] <- "Variable"
  
  for(j in 1:nrow(genes)){
    temp = sub3
    
    temp$genotype = ifelse(temp$Gene  == genes[j,1], "Mut", "WT")
    
    temp = unique(select(temp, Sample, Variable, genotype))
    
    n = as.numeric(length(unique(temp$Sample)))
    
    n_check = as.numeric(length(which(temp$genotype == "Mut")))
    
    if(n_check >= 5){
      print(j)
      # calculate a p-value and effect size for the difference in clinical feaures based on genotype differences
      p_val =  wilcox.test(as.numeric(temp$Variable) ~ temp$genotype, alternative = "two.sided")$p.value
      
      temp$Variable = as.numeric(temp$Variable)
      
      effect_size = as.numeric(cohens_d(Variable ~ genotype, data = temp)$effsize)
      n_n1 = as.numeric(length(which(temp$genotype == "Mut")))
      n_n2 = as.numeric(length(which(temp$genotype == "WT")))
      effect_size_ci = cohen.d.ci(d = effect_size, n = n, n1 = n_n1, n2 = n_n2)
      
      
      variable_stats <- data.frame(matrix(NA, nrow = 1, ncol = 7))
      names(variable_stats) <- c("Variable", "Gene", "effect_size", "CI_lower", "CI_upper", "p_value", "n_mut")
      
      variable_stats[1,1] <- paste(var)
      variable_stats[1,2] <- as.character(genes[j,1])
      variable_stats[1,3] <- effect_size
      variable_stats[1,4] <- effect_size_ci[1]
      variable_stats[1,5] <- effect_size_ci[3]
      variable_stats[1,6] <- p_val
      variable_stats[1,7] <- n_check
      
      
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


# gene_pairs_list_final_adj$point_label = paste(gene_pairs_list_final_adj$Gene_1, " + ", gene_pairs_list_final_adj$Gene_2, sep = "")
gene_pairs_list_final_adj$point_label = gene_pairs_list_final_adj$Gene


for(i in 1:nrow(gene_pairs_list_final_adj)){
  if(gene_pairs_list_final_adj$q_val[i] > 0.05){
    gene_pairs_list_final_adj$point_label[i] = ""
  }
}

gene_pairs_list_final_adj$sig_color = as.factor(gene_pairs_list_final_adj$sig_color)

max_p = max(-log10(gene_pairs_list_final_adj$p_value))

gene_pairs_list_final_adj$Variable = factor(gene_pairs_list_final_adj$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM_blast_percent', 'PB_blast_percent', 'Age'))

# plot results per variable
p = ggplot(gene_pairs_list_final_adj, aes(x=gene_pairs_list_final_adj$effect_size, y=-log10(gene_pairs_list_final_adj$p_value), color = sig_color, size = n_mut)) +
  # geom_label_repel(aes(label=point_label),box.padding = .5, hjust=0, vjust=1.5, size = 2.5) +
  geom_hline(yintercept = 2.430129,  linetype = "dashed", color = "lightgrey") +
  geom_point(alpha = 0.75) +
  geom_point(shape = 21, color = "black", alpha = 0.25)  +
  scale_size_area(max_size = 10,breaks=c(10,25,50,100,250,350)) +
  geom_label_repel(aes(label=point_label), force = 25, size = 2.5) +
  # geom_hline(yintercept = 3.563053,  linetype = "dashed", color = "lightgrey") +
  theme_cowplot() +
  scale_colour_manual(values = c("sig_pos_ES"="#b35806", "sig_neg_ES" = "#542788", "not_sig_ES"="lightgrey")) +
  theme(legend.position="top") +
  # geom_hline(yintercept = 1.30103,  linetype = "dashed", color = "darkgrey") +
  ylab(label= "-log10(p-value)") +
  xlab(label= "Effect Size (Mut vs. WT)") +
  # scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(color="black", size=20)) 

g = guide_legend(override.aes=list(colour="lightgrey"), "n. mut")

p + facet_wrap(. ~ Variable, ncol = 1, nrow = 7) + guides(size = g, color = FALSE) +
  theme(
    # legend.title = element_blank(), 
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid")) 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Individual_genotype_clinical_correlations.png", dpi = 300, width = 3.75, height = 15, units = "in")

write.csv(gene_pairs_list_final_adj, file = "~/Desktop/MetaAML_results/Figure_2/Supplimental/gene_clinical_features.csv")


 # forrest plot ####
# individual mutation's HRs

temp_final_hr_2$gene <- factor(temp_final_hr_2$gene, levels = temp_final_hr_2$gene[order(temp_final_hr_2$HR)])
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
  if(temp_final_hr_2$gene[i] %in% DNA_methylation){
    temp_final_hr_2$mutation_category[i] <- "DNA Methylation"
  }
  if(temp_final_hr_2$gene[i] %in% Chromatin_cohesin){
    temp_final_hr_2$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(temp_final_hr_2$gene[i] %in% RTK_RAS_Signaling){
    temp_final_hr_2$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(temp_final_hr_2$gene[i] %in% Splicing){
    temp_final_hr_2$mutation_category[i] <- "Splicing"
  }
  if(temp_final_hr_2$gene[i] %in% Transcription){
    temp_final_hr_2$mutation_category[i] <- "Transcription"
  }
  if(temp_final_hr_2$gene[i] == "NPM1"){
    temp_final_hr_2$mutation_category[i] <- "NPM1"
  }
  if(temp_final_hr_2$gene[i] %in% Tumor_suppressors){
    temp_final_hr_2$mutation_category[i] <- "Tumor suppressors"
  }
}

# extract the p-value and hazard ratio for the individual interactions
p = round(as.numeric(summary(model)$sctest[3]), 3)

p = ifelse(p < 0.001, paste0("p < 0.001"), paste("p =", p))

hr = paste("HR = ", round(as.numeric(forest$HR), 2), " (", round(as.numeric(forest$lower_95), 2), "-", round(as.numeric(forest$upper_95), 2), ")", sep = "")

p_hr = paste(p, "; ", hr, sep = "")


# color coded by mutation category plot
ggplot(temp_final_hr_2, aes(x = reorder(gene, -HR), y = HR, label = temp_final_hr_2$p_text)) +
  geom_hline(yintercept=.5, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=3, linetype="dashed", color = "#d9d9d9") +
  geom_text(aes(gene, upper_95), hjust = 0, nudge_y = 0.35, size = 3) +
  theme_cowplot() +
  ylim(0,5.25) +
  geom_pointrange(size = .75, stat = "identity", shape = 15, 
                  aes(x = gene, ymin = lower_95, ymax = upper_95, y = HR, color = mutation_category)) +
  scale_color_manual(name = "", values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF")) +
  ylab("Hazard Ratio")+
  xlab(NULL) +
  theme(legend.position = c(0.6,0.25),
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/individual_gene_HR_forrest_plot_de_novo.pdf", dpi = 300, width = 9, height = 6, units = "in")

