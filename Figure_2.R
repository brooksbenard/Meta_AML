# ========================================================================================================================================== #
# Figure_2.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 08/23/2021
# Description: This script will perform survival and mutation co-occurence analyses as seen in Figure 2 and related suppliments of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# ================ #
# Load packages ####
# ================ #
# Package names
packages <- c("ggplot2", "dplyr", "cometExactTest", "stringr", "maditr", "reshape2", "data.table", "epitools", "corrplot", "plyr", "muhaz", "reshape", "survival", "survivalAnalysis", "survMisc", "survminer", "ggsci", "vegan", "ggrepel", "ggforce", "rstatix", "effsize", "psych", "devtools", "ggrepel", "cowplot", "RColorBrewer", "readr")

# devtools::install_github("thomasp85/ggforce")
library("ggforce")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# create directories
dir.create("~/Desktop/MetaAML_results/Figure_2")
dir.create("~/Desktop/MetaAML_results/Figure_2/Supplimental")

# load from data frame created in the Figure_1 script
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

#### Individual and Pairwise genotype associations with clinical features ####
sub = subset(final_data_matrix, mut_freq_gene > 50 & final_data_matrix$Subset == "de_novo")

# make sure that the FLT3 symbols are annotated well
for(i in 1:nrow(sub)){
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] %in% c("ITD", "INDEL")){
    sub$Gene[i] <- "FLT3-ITD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] %in% c("SNV", "Deletion", "other")){
    sub$Gene[i] <- "FLT3-TKD"
  }
}

genes = data.frame(unique(sub$Gene))

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
age2 = ggplot(sub, aes(x=reorder(Gene, -Age, median, na.rm = T), y=Age)) +
          geom_boxplot(fill="#bdbdbd", color="black") +
          theme_cowplot(font_size = 15) +
          theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                           size = 10, hjust = .5)) +
          xlab(label = NULL) +
          ylab("Age") +
          geom_hline(yintercept=mean(sub$Age, na.rm = T), linetype="dashed", color = "#99000d")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Age_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

# Platelet
platelet_2 = ggplot(sub, aes(x=reorder(Gene, -Platelet, median, na.rm = T), y=Platelet)) +
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
heme_2 = ggplot(sub, aes(x=reorder(Gene, -Hemoglobin, median, na.rm = T), y=Hemoglobin)) +
              geom_boxplot(fill="#fe9929", color="black") +
              theme_cowplot(font_size = 15) +
              theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                               size = 10, hjust = .5)) +
              xlab(label = NULL) +
              ylab("Hemoglobin (g/dL)") +
              geom_hline(yintercept=mean(sub$Hemoglobin, na.rm = T), linetype="dashed", color = "#99000d") 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/Hemoglobin_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

# WBC
wbc2 = ggplot(sub, aes(x=reorder(Gene, -WBC, median, na.rm = T), y=WBC)) +
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
ldh2 = ggplot(sub, aes(x=reorder(Gene, -LDH, median, na.rm = T), y=LDH)) +
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
bm_blast = ggplot(sub, aes(x=reorder(Gene, -BM_blast_percent, median, na.rm = T), y=BM_blast_percent)) +
                geom_boxplot(fill="#e6ab02", color="black") +
                theme_cowplot(font_size = 15) +
                theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                                 size = 10, hjust = .5)) +
                xlab(label = NULL) +
                ylab("BM Blast %") +
                geom_hline(yintercept=mean(sub$BM_blast_percent, na.rm = T), linetype="dashed", color = "#99000d") 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/BM_blast_percent_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 

# PB_blast_percent
pb_blast = ggplot(sub, aes(x=reorder(Gene, -PB_blast_percent, median, na.rm = T), y=PB_blast_percent)) +
                  geom_boxplot(fill="#8dd3c7", color="black") +
                  theme_cowplot(font_size = 15) +
                  theme(axis.text.x = element_text(angle = 45, vjust = .75, 
                                                   size = 10, hjust = .5)) +
                  xlab(label = NULL) +
                  ylab("PB Blast %") +
                  geom_hline(yintercept=mean(sub$PB_blast_percent, na.rm = T), linetype="dashed", color = "#99000d")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/PB_blast_percent_vs_gene.png", dpi = 300, width = 10, height = 2.5, units = "in") 


# put all of the plots into a panel
ggarrange(wbc2, heme_2, platelet_2, ldh2, bm_blast, pb_blast, age2,
          nrow = 7,
          ncol = 1,
          align = "hv"
)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/all_variables_panel.pdf", dpi = 300, width = 10, height = 17.5, units = "in")




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
      effect_size_ci = psych::cohen.d.ci(d = effect_size, n = n, n1 = n_n1, n2 = n_n2)
      
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
  # geom_point(shape = 21, color = "black", alpha = 0.25)  +
  geom_label_repel(aes(label=point_label), size = 2, force = 75, max.overlaps = 5) +
  scale_size_area(max_size = 10,breaks=c(10,25,50,100,250,350)) +
  scale_colour_manual(values = c("sig_pos_ES"="#b35806", "sig_neg_ES" = "#542788", "not_sig_ES"="lightgrey")) +
  theme(legend.position="right") +
  ylab(label= "-log10(p-value)") +
  xlab(label= "Effect Size (co-mut vs. others)") +
  theme(plot.title = element_text(color="black", size=20)) 

g = guide_legend(override.aes=list(size="lightgrey"), title = "n. co-mut")

p = p + facet_wrap(. ~ Variable, ncol = 7) + ggplot2::guides(color = "none", size = g) + 
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid"))

ggsave(plot = p, filename = "~/Desktop/MetaAML_results/Figure_2/Pairwise_genotype_clinical_correlations.pdf", dpi = 300, width = 15, height = 5, units = "in")



#### statistical analysis of mutation co-occurence ####
# subset analysis to the de novo cohort
final_data_matrix_2_sub <- subset(final_data_matrix, Subset == "de_novo")

# make sure that the FLT3 symbols are the same
for(i in 1:nrow(final_data_matrix_2_sub)){
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] %in% c("ITD", "Deletion")){
    final_data_matrix_2_sub$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] %in% c("SNV", "INDEL", "other")){
    final_data_matrix_2_sub$Gene[i] <- "FLT3-TKD"
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
# select informative columns
final_data_matrix_2_sub = select(final_data_matrix_2_sub, Sample, Gene, VAF)

# transform the data frame
final_data_matrix_2_sub <- dcast(final_data_matrix_2_sub, Sample ~ Gene, value.var="VAF")

rownames(final_data_matrix_2_sub) <- final_data_matrix_2_sub$Sample
final_data_matrix_2_sub$Sample <- NULL
final_data_matrix_2_sub[final_data_matrix_2_sub != 0] <- 1
temp_for_odds_ratio = final_data_matrix_2_sub
final_data_matrix_2_sub <- as.data.frame(t(final_data_matrix_2_sub))

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

p1 = ggplot(df, aes(x=log(odds_ratio), y=-log10(fishers_q), color = significant, size = n_cooccur)) +
  theme_cowplot() +
  geom_hline(yintercept=1.30103, linetype="dashed", color = "#d9d9d9") +
  geom_point(alpha = .75) +
  geom_label_repel(aes(label=labels),nudge_y = 5, size = 3, force = 25, max.overlaps = 5) +
  scale_colour_manual(values = c("Co-occuring"= "#b2182b", "Mutually exclusive"="#2166ac", "NS"="lightgrey")) + 
  # theme(legend.position = c(0.05,.85)) +
  ylab(label= "-log10(q-value)") +
  xlab(label= "log(Odds Ratio)") +
  scale_size_area(max_size = 10,breaks=c(25,50,100,200,300)) +
  labs(title = NULL) 

g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p1 = p1 + guides(size = "none", color = "none")  +
  annotate(geom="text", x = 1.9, y = 0, label="q < 0.05",
           color="grey")

 ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/volcano_plot_co_occurence.pdf", dpi = 300, width = 5, height =  5, units = "in")


# co-occurence and survival ####

# load data and subset to desired cohort
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_sub = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_sub = subset(final_data_matrix_sub, final_data_matrix_sub$mut_freq_gene >= 100)
final_data_matrix_sub$Time_to_OS = (final_data_matrix_sub$Time_to_OS/365)

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(final_data_matrix_sub)){
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] %in% c("ITD", "INDEL")){
    final_data_matrix_sub$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_sub$Gene[i] == "FLT3" & final_data_matrix_sub$variant_type[i] %in% c("SNV", "Deletion", "other")){
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

p2 = ggplot(df, aes(x=df$HR, y=-log10(df$q_value), color = factor(significant), size = num_pts)) +
  theme_cowplot() +
  geom_hline(yintercept=1, linetype="dashed", color = "#d9d9d9") +
  geom_point(alpha = .75) +
  scale_size_area(max_size = 10,breaks=c(25,50,100,200,300)) +
  geom_label_repel(aes(label=labels),size = 3, force = 25, max.overlaps = 1) +
  scale_colour_manual(values = c("Favorable"= "#1b7837", "Unfavorable" = "#762a83", "NS"="grey")) + 
  theme_cowplot() +
  # theme(legend.position = c(0.05,.85)) +
  ylab(label= "-log10(q-value)") +
  xlab(label= "Hazard Ratio") +
  labs(title = NULL) +
  theme(plot.title = element_text(color="black", size=20)) 
g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p2 = p2  + guides(size = g, color = "none")  +
  annotate(geom="text", x = 1.25, y = .75, label="q < 0.05",
           color="grey")+
  theme(legend.position = "none") 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/volcano_plot_co_occurence_survival.pdf", dpi = 300, width = 5, height =  5, units = "in")



# scattterplot of HR and odds ratio ####
# read in the results from MetaAML_co_occurence_mutual_exlclusive_analysis.R
odds_ratio = read.csv("~/Desktop/MetaAML_results/Figure_2/Fishers/odds_ratio_and_fishers_results.csv")
odds_ratio = select(odds_ratio, gene1, gene2, odds_ratio, fishers_q)
colnames(odds_ratio) = c("gene_1", "gene_2", "odds_ratio", "q_value_odds")

hazard_ratio = read.csv("~/Desktop/MetaAML_results/Data/Tables/pairwise_mutations_hazard_ratio.csv")
hazard_ratio = hazard_ratio %>%
  select(gene_1, gene_2, HR, q_value, num_pts)
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

p3 = ggplot(hr_odds, aes(x = log(odds_ratio), y = log(HR), color = as.factor(color), size = num_pts)) +
  geom_point(alpha = 0.75) +
  xlab("log(Odds Ratio)") +
  ylab("log(Hazard Ratio)") +
  scale_color_manual(values = c("#1b7837",  "#762a83", "grey")) +
  scale_fill_manual(values = c("#1b7837",  "#762a83", "grey")) +
  geom_hline(yintercept=log(1), linetype="dashed", color = "grey") +
  geom_vline(xintercept=log(1), linetype="dashed", color = "grey") +
  geom_smooth(method = lm, aes(fill = color),  alpha = .25) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -1.95, 
    label.y = c(-.14,-.16,-.18)
    
  ) +
  theme_cowplot() +
  scale_size_area(max_size = 10,breaks=c(25,50,100,200,300)) +
  theme(plot.title = element_text(color="black", size=20))
g = guide_legend(override.aes=list(colour="lightgrey"), "n. co-mut")

p3  = p3 + guides(size = "none", color = "none", fill = "none") +
  annotate(geom="text", x = 1.8, y = .25, label="26%",
           color="#762a83") +
  annotate(geom="text", x = -1.8, y = .25, label="28%",
           color="#762a83")+
  annotate(geom="text", x = 1.8, y = -.2, label="12%",
           color="#1b7837") +
  annotate(geom="text", x = -1.8, y = -.2, label="14%",
           color="#1b7837")


 ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/correlation_of_odds_ratio_and_HR_de_novo_2.pdf", dpi = 300, width = 5, height = 5, units = "in")

  write.csv(hr_odds, "~/Desktop/MetaAML_results/Data/Tables/correlation_of_odds_ratio_and_HR_de_novo.csv")

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



# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
p = ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="right", vjust = -0.8)

g = guide_legend(override.aes=list(colour="lightgrey"), "n. nco-mut")

p  + guides(size = g, color = "none", fill = "none") 
      
ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/correlation_of_odds_ratio_and_HR_de_novo_all.pdf", dpi = 300, width = 15, height = 5, units = "in")


# tripple mutation analysis ####
# analyze the distribution and association of outcomes with triple mutation subsets

load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# count the frequencies of recurrent mutations and select only genes with n > 50
final_data_matrix_2_alt = final_data_matrix %>%
  subset(Subset == "de_novo")

for(i in 1:nrow(final_data_matrix_2_alt)){
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "ITD"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "SNV"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "Deletion"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "INDEL"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "other"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-TKD"
  }
}

three_mutations = final_data_matrix_2_alt %>%
  group_by(Gene) %>%
  dplyr::count(Gene, sort = TRUE) %>%
  filter(n > 50)

# find all 3-way mutation combinations
genes = as.data.frame(t(combn(as.vector(unique(three_mutations$Gene)),3)))
genes$n_3_mut = NA

# select columns for analysis from data matrix
three_mutations_clinical = final_data_matrix_2_alt %>%
  select(Sample, Gene, Time_to_OS, Censor)

# select individual patients
three_mutations_pts = final_data_matrix_2_alt %>%
  select(Sample) %>%
  distinct() %>%
  mutate(three_mutations = 0)

mut_list = list()

for(i in 1:nrow(genes)){
  print(i)
  gene_list = c(genes[i,1:3])
  for(j in 1:nrow(three_mutations_pts)){
    pt = three_mutations_clinical %>%
      subset(Sample == three_mutations_pts$Sample[j])
    if(gene_list$V1[1] %in% pt$Gene & gene_list$V2[1] %in% pt$Gene & gene_list$V3[1] %in% pt$Gene){
      # print("TRUE")
      three_mutations_pts$three_mutations[j] = 1
    } else {
      # print("FALSE")
      three_mutations_pts$three_mutations[j] = 0
    }
  }
  genes$n_3_mut[i] = length(which(three_mutations_pts$three_mutations == 1))
}

genes_all_3_mut_n_over_5 = genes %>%
  filter(n_3_mut >= 10)

# now perform survival analysis between genotypes in the set of mutations with at least 10 patiens with 3 co-occuring mutations

# create columns for the 3 gene genotypes
genes_all_3_mut_n_over_5$genotype = do.call(paste, c(genes_all_3_mut_n_over_5[,1:3], sep = '_'))

cnames = c(dcast(genes_all_3_mut_n_over_5,.~genotype))
cnames[1] = NULL

for(i in cnames){
  three_mutations_pts[,i] = 0
}

three_mutations_pts <- three_mutations_pts %>%
  mutate_if(is.logical, as.numeric)

three_mutations_pts[three_mutations_pts == 0] <- "WT"

# function to find genes not in a list
"%ni%" <- Negate("%in%") 

# loop through all of the mutation genotype combinations and label each patient based on the genotype. This will be done foe specific genotypes and broad relationships
three_mutations_pts_simple = three_mutations_pts

for(i in 1:nrow(genes_all_3_mut_n_over_5)){
  print(i)
  gene_list = c(genes_all_3_mut_n_over_5[i,1:3])
  mut_genotype = genes_all_3_mut_n_over_5[i,5]
  c_num = as.numeric(which(colnames(three_mutations_pts) == mut_genotype))
  for(j in 1:nrow(three_mutations_pts)){
    pt = three_mutations_clinical %>%
      subset(Sample == three_mutations_pts$Sample[j])
    if(gene_list$V1[1] %in% pt$Gene & gene_list$V2[1] %in% pt$Gene & gene_list$V3[1] %in% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V1[1], gene_list$V2[1], gene_list$V3[1], sep = "_")
      three_mutations_pts_simple[j,c_num] = "3_mut"
    } 
    if(gene_list$V1[1] %in% pt$Gene & gene_list$V2[1] %in% pt$Gene & gene_list$V3[1] %ni% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V1[1], gene_list$V2[1], sep = "_")
      three_mutations_pts_simple[j,c_num] = "2_mut"
    }
    if(gene_list$V1[1] %in% pt$Gene & gene_list$V2[1] %ni% pt$Gene & gene_list$V3[1] %in% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V1[1], gene_list$V3[1], sep = "_")
      three_mutations_pts_simple[j,c_num] = "2_mut"
    }
    if(gene_list$V1[1] %ni% pt$Gene & gene_list$V2[1] %in% pt$Gene & gene_list$V3[1] %in% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V2[1], gene_list$V3[1], sep = "_")
      three_mutations_pts_simple[j,c_num] = "2_mut"
    }
    if(gene_list$V1[1] %in% pt$Gene & gene_list$V2[1] %ni% pt$Gene & gene_list$V3[1] %ni% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V1[1])
      three_mutations_pts_simple[j,c_num] = "1_mut"
    }
    if(gene_list$V1[1] %ni% pt$Gene & gene_list$V2[1] %in% pt$Gene & gene_list$V3[1] %ni% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V2[1])
      three_mutations_pts_simple[j,c_num] = "1_mut"
    }
    if(gene_list$V1[1] %ni% pt$Gene & gene_list$V2[1] %ni% pt$Gene & gene_list$V3[1] %in% pt$Gene){
      three_mutations_pts[j,c_num] = paste(gene_list$V3[1])
      three_mutations_pts_simple[j,c_num] = "1_mut"
    }
  }
}

# select the survival data
three_mutations_clinical_surv = three_mutations_clinical %>%
  select(Sample, Time_to_OS, Censor) %>%
  unique()

# pair survival data with genotypes on a granular level
genotype_survival = left_join(three_mutations_pts, three_mutations_clinical_surv, by = "Sample")
genotype_survival$Censor = as.numeric(genotype_survival$Censor)
genotype_survival$Time_to_OS = genotype_survival$Time_to_OS/365

write_csv(genotype_survival, "~/Desktop/MetaAML_results/Figure_2/genotype_combinations.csv")

# pair survival data with genotypes on a broad level
genotype_survival_simple = left_join(three_mutations_pts_simple, three_mutations_clinical_surv, by = "Sample")
genotype_survival_simple$Censor = as.numeric(genotype_survival_simple$Censor)
genotype_survival_simple$Time_to_OS = genotype_survival_simple$Time_to_OS/365

hr_info = list()
n = 1
options(scipen = 999)
for(i in 3:50){
  
  genotype_survival_simple[[i]][!grepl("2_mut|3_mut",genotype_survival_simple[[i]])]<-NA
  
  genotype_names1 = paste(colnames(genotype_survival_simple)[i])
  genotype_names = gsub("_", " + ", genotype_names1, fixed=TRUE)
  
  fit <- survfit(Surv(Time_to_OS, Censor) ~ genotype_survival_simple[[i]], data = genotype_survival_simple)
  names(fit$strata) = as.character(names(fit$strata))
  names(fit$strata) = gsub("genotype_survival_simple[[i]]=", "", names(fit$strata), fixed = T)
  names(fit$strata) = gsub("_", " ", names(fit$strata), fixed = T)
  
  surv_plot = ggsurvplot(fit,
                         conf.int = F,
                         risk.table = T,
                         pval = T,
                         pval.coord = c(0, 0),
                         title = paste(genotype_names),
                         palette = c("darkgrey", "darkred"),
                         risk.table.y.text = FALSE,
                         tables.y.text = FALSE,
                         legend = "right",
                         legend.title = "",
                         xlab = "Years")
  
  print(surv_plot)
  png(filename = paste("~/Desktop/MetaAML_results/Figure_2/Supplimental/2_3_genoptypes/broad/",genotype_names1,"survival.png", sep = ""), res = 300, width = 5, height = 5, units = "in")
  # ggsave(file = paste("~/Desktop/Majeti_Lab/Manuscripts/Meta_AML/Final_submission_materials/Nat_Comm/Revisions/Comment_2/2_3_genotypes/broad/",genotype_names,"survival.pdf", sep = ""), print(surv_plot))
  
  surv_plot
  print(surv_plot)
  dev.off()
  
  
  # extract data from survfit
  model <- coxph( Surv(Time_to_OS, Censor) ~ genotype_survival_simple[[i]],
                  data = genotype_survival_simple )
  
  # extract the informative data from the survival model
  array_dat = summary(model)$conf.int[1:4]
  array_dat[5] = genotype_names
  
  # extract the log-rank p-value for the individual comparisons
  array_dat[6] = summary(model)$sctest[3]
  array_dat = array_dat[-2]
  
  forest_plot_data <- data.frame("Genes" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = round(as.numeric(array_dat[2]),2), "Upper_CI" = round(as.numeric(array_dat[3]),2), "log_rank_p" = format(array_dat[5], scientific=F))
  
  hr_info[[n]] <- forest_plot_data
  n=n+1   
  
}
hr_info_final = as.data.frame(do.call(rbind, hr_info)) 
hr_info_final = unique(hr_info_final)
hr_info_final$Genes <- factor(hr_info_final$Genes, levels = hr_info_final$Genes[order(hr_info_final$HR)])

write_csv(hr_info_final, "~/Desktop/MetaAML_results/Figure_2/3_mut_hr_table.csv")


hr_info_final$sig_color = 0

for(i in 1:nrow(hr_info_final)){
  if(hr_info_final$log_rank_p[i] < 0.05 & hr_info_final$HR[i] < 1){
    hr_info_final$sig_color[i] = 1
  }
  if(hr_info_final$log_rank_p[i] < 0.05 & hr_info_final$HR[i] > 1){
    hr_info_final$sig_color[i] = 2
  }
}

hr_info_final$sig_color = as.factor(hr_info_final$sig_color)
hr_info_final$fdr = p.adjust(hr_info_final$log_rank_p, method = "fdr")
hr_info_final$sig_color = as.factor(hr_info_final$sig_color)
hr_info_final = subset(hr_info_final, hr_info_final$Upper_CI != "Inf")
hr_info_final$categories <- reorder(hr_info_final$categories, hr_info_final$HR)


hr_info_final$p_text = NA
for(i in 1:nrow(hr_info_final)){
  if(hr_info_final$log_rank_p[i] < 0.05){
    hr_info_final$p_text[i] = hr_info_final$log_rank_p[i]
  }
}
hr_info_final$q_text = NA
for(i in 1:nrow(hr_info_final)){
  if(hr_info_final$fdr[i] < 0.5){
    hr_info_final$q_text[i] = hr_info_final$fdr[i]
  }
}

hr_info_final$p_text = as.numeric(hr_info_final$p_text)
hr_info_final$q_text = as.numeric(hr_info_final$q_text)

hr_info_final$p_text = round(hr_info_final$p_text, 3)
hr_info_final$q_text = round(hr_info_final$q_text, 2)

for(i in 1:nrow(hr_info_final)){
  if(hr_info_final$log_rank_p[i] < 0.01){
    hr_info_final$p_text[i] = "p < 0.01"
  }
  if(hr_info_final$log_rank_p[i] >= 0.01 & hr_info_final$log_rank_p[i] <= 0.05){
    hr_info_final$p_text[i] = paste("p =", paste(hr_info_final$p_text[i]))
  }
  if(hr_info_final$fdr[i] < 0.01){
    hr_info_final$q_text[i] = "q < 0.01"
  }
  if(hr_info_final$fdr[i] >= 0.01 & hr_info_final$log_rank_p[i] <= 0.05){
    hr_info_final$q_text[i] = paste("q =", paste(hr_info_final$q_text[i]))
  }
}

hr_info_final$p_q_text = paste(hr_info_final$p_text, "; ", hr_info_final$q_text, sep = "")

for(i in 1:nrow(hr_info_final)){
  if(hr_info_final$log_rank_p[i] > 0.05){
    hr_info_final$p_q_text[i] = ""
  }
  if(hr_info_final$log_rank_p[i] > 0.05){
    hr_info_final$p_q_text[i] = ""
  }
}

hr_info_final$HR = as.numeric(hr_info_final$HR)
hr_info_final$Lower_CI = as.numeric(hr_info_final$Lower_CI)
hr_info_final$Upper_CI = as.numeric(hr_info_final$Upper_CI)

hr_info_final$Genes <- factor(hr_info_final$Genes, levels = hr_info_final$Gene[order(hr_info_final$HR)])

names(hr_info_final)[1] = "genotype"

genes_all_3_mut_n_over_5$genotype =gsub("_", " + ", genes_all_3_mut_n_over_5$genotype, fixed=TRUE)
hr_info_final$genotype = gsub("_", " + ", hr_info_final$genotype, fixed=TRUE)

hr_info_final = left_join(hr_info_final, genes_all_3_mut_n_over_5, by = "genotype")

hr_info_final$genotype <- reorder(hr_info_final$genotype, hr_info_final$HR)

ggplot(hr_info_final, aes(x = genotype, y = HR, label = p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(genotype, Upper_CI), hjust = 0, nudge_y = 0.6) +
  geom_point(aes(size = n_3_mut, color = sig_color), shape = 19, alpha = .9) +
  geom_segment(data = hr_info_final, aes(y = Lower_CI, yend = Upper_CI, x = genotype, xend = genotype, col= sig_color), size = 0.5) +
  scale_color_manual(values = c("1" = "#1b7837", "2" = "#762a83"))+
  ylab("Hazard Ratio\n(3 muts vs. 2 muts)")+
  ylim(0,8)+
  theme_cowplot() +
  theme(legend.position = c(0.7, .1),
        axis.title.y=element_blank()) +
  coord_flip() +
  guides(color = FALSE,
         size =  guide_legend(override.aes=list(colour="lightgrey"), title = "# of patients with\n3 mutations"))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/_3_mut_forrest_plot.pdf", dpi = 300, width = 8.5, height = 10, units = "in")


# distribution of genotype frequencies
# genotypes with enough survival data for analysis
genotype_list = hr_info_final %>%
  subset(sig_color != 0) %>%
  select(genotype)

genes_all_3_mut_n_over_5 = genes_all_3_mut_n_over_5 %>%
  mutate(`Survival association` = ifelse(genotype %in% genotype_list$genotype, "Yes", "No")) 


ggplot(data=genes_all_3_mut_n_over_5, aes(x=reorder(genotype, -n_3_mut), y = n_3_mut, fill = `Survival association`)) +
  theme_cowplot(font_size = 15) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('grey', "darkred")) +
  geom_text(aes(label=n_3_mut), vjust=-0.3, size=5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = c(.9, 0.5)) +
  ylab("\n\n\n# of patiens\nwith 3 mutations") +
  ylim(0,140)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/three_mut_genotype_distribution.pdf", dpi = 300, width = 20, height = 5, units = "in")



# Supplimental ####
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


# multi-mut per gene/gategory ####
# 4. The authors should make it more clear how they are addressing cases in which there are multiple mutations in the same gene or pathway (i.e. RAS/MAPK-pathway mutations). Does convergent evolution influence clinical correlates or survival outcomes?
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

final_data_matrix_2_alt = final_data_matrix %>%
  subset(Subset == "de_novo")

for(i in 1:nrow(final_data_matrix_2_alt)){
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "ITD"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "SNV"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "Deletion"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "INDEL"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2_alt$Gene[i] == "FLT3" & final_data_matrix_2_alt$variant_type[i] == "other"){
    final_data_matrix_2_alt$Gene[i] <- "FLT3-TKD"
  }
}

# add a column for annotating the mutation category
final_data_matrix_2_alt$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(final_data_matrix_2_alt)){
  if(final_data_matrix_2_alt$Gene[i] %in% DNA_methylation){
    final_data_matrix_2_alt$mutation_category[i] <- "DNA Methylation"
  }
  if(final_data_matrix_2_alt$Gene[i] %in% Chromatin_cohesin){
    final_data_matrix_2_alt$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(final_data_matrix_2_alt$Gene[i] %in% RTK_RAS_Signaling){
    final_data_matrix_2_alt$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(final_data_matrix_2_alt$Gene[i] %in% Splicing){
    final_data_matrix_2_alt$mutation_category[i] <- "Splicing"
  }
  if(final_data_matrix_2_alt$Gene[i] %in% Transcription){
    final_data_matrix_2_alt$mutation_category[i] <- "Transcription"
  }
  if(final_data_matrix_2_alt$Gene[i] == "NPM1"){
    final_data_matrix_2_alt$mutation_category[i] <- "NPM1"
  }
  if(final_data_matrix_2_alt$Gene[i] %in% Tumor_suppressors){
    final_data_matrix_2_alt$mutation_category[i] <- "Tumor suppressors"
  }
}

convergent_gene = final_data_matrix_2_alt %>%
  group_by(Gene) %>%
  dplyr::count(Sample, sort = TRUE) 

convergent_clinical = final_data_matrix_2_alt %>%
  subset(Subset == "de_novo") %>%
  select(Sample, BM_blast_percent, PB_blast_percent, WBC, Hemoglobin, LDH, Platelet, Age, Censor, Time_to_OS) %>%
  unique()

convergent_clinical_merged = left_join(convergent_clinical, convergent_gene, by = "Sample")

convergent_genes = final_data_matrix_2_alt %>%
  group_by(Gene) %>%
  dplyr::count(Sample, sort = TRUE) %>%
  filter(n > 1) %>%
  distinct(Gene)

convergent_clinical_merged <- convergent_clinical_merged %>% 
  mutate_at(c(2:10), as.numeric)

variable_df_list_1 = list()
variable_df_list_2 =  list()
forest_plot_data_list = list()

k = 1
for(i in 1:nrow(convergent_genes)){
  sub_dat_final = convergent_clinical_merged %>%
    subset(Gene == convergent_genes$Gene[i]) %>%
    distinct(Sample, BM_blast_percent, PB_blast_percent, WBC, Hemoglobin, LDH, Platelet, Age, Censor, Time_to_OS, n, Gene)
  
  sub_dat_final$n = as.factor(ifelse(sub_dat_final$n>1, "2", "1"))
  
  if(length(which(sub_dat_final$n != 1)) > 9){
    l = 1
    for(j in 2:8){
      # calculate a p-value and effect size for the difference in clinical feaures based on number of mutations
      # need to dynamically chage the column name for the variable of interest in order to calculate the effect size. will change back
      col_name = as.character(colnames(sub_dat_final[j])[j])
      names(sub_dat_final)[names(sub_dat_final) == col_name] <- "Variable"
      
      sub_dat_final[[j]] = as.numeric(sub_dat_final[[j]])
      
      # calculate effect sizes for clinical correlates
      p_val =  wilcox.test(sub_dat_final[[j]] ~ sub_dat_final$n, alternative = "two.sided")$p.value
      
      effect_size = cohens_d(data = sub_dat_final, Variable ~ n)$effsize
      effect_size_ci = cohens_d(data = sub_dat_final, Variable ~ n, ci = TRUE)
      n_n1 = as.numeric(length(which(sub_dat_final$n == 1)))
      n_n2 = as.numeric(length(which(sub_dat_final$n == 2)))
      effect_size_ci = psych::cohen.d.ci(d = effect_size, n = n, n1 = n_n1, n2 = n_n2)
      
      variable_df <- data.frame(matrix(NA, nrow = 1, ncol = 6))
      names(variable_df) <- c("Variable", "Mutated_Gene", "effect_size", "CI_lower", "CI_upper", "p_value")
      
      variable_df[1,1] <- paste(col_name)
      variable_df[1,2] <- convergent_genes$Gene[i]
      variable_df[1,3] <- effect_size
      variable_df[1,4] <- effect_size_ci[1,1]
      variable_df[1,5] <- effect_size_ci[1,3]
      variable_df[1,6] <- p_val
      
      # Add each list in the loop to a list of lists and rename column
      variable_df_list_1[[l]] = variable_df 
      l = l + 1 
      
      names(sub_dat_final)[names(sub_dat_final) == "Variable"] = col_name
      
    }
    variable_df_list_1_temp = do.call(rbind, variable_df_list_1)
    variable_df_list_1_temp = variable_df_list_1_temp %>%
      mutate(q_val = p.adjust(variable_df_list_1_temp$p_value))
    
    variable_df_list_2[[k]] = variable_df_list_1_temp
    k=k+1
    
    # calculate hazard ratios and confidence intervals for survival outcomes
    sub_dat_final$Censor = as.numeric(sub_dat_final$Censor)
    # run the Cox model
    model <- coxph( Surv(Time_to_OS, Censor) ~ n,
                    data = sub_dat_final )

    # extract the informative data from the survival model
    array_dat = summary(model)$conf.int[1:4]
    array_dat[5] = convergent_genes$Gene[i]
    
    # extract the log-rank p-value for the individual comparisons
    array_dat[6] = summary(model)$sctest[3]
    array_dat = array_dat[-2]
    
    forest_plot_data <- data.frame("Variable" = "Survival", "Mutated_Gene" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = round(as.numeric(array_dat[2]),2), "Upper_CI" = round(as.numeric(array_dat[3]),2), "log_rank_p" = array_dat[5])
    
    forest_plot_data_list[[i]] = forest_plot_data
    
  }
  
}
# bind all of the clinical outcomes data and correct for FDR
forest_plot_data_list_all = do.call(rbind, forest_plot_data_list)
forest_plot_data_list_all$q_val = p.adjust(forest_plot_data_list_all$log_rank_p)

# bind all of the clinical correlates data
variable_df_list_all <- do.call(rbind, variable_df_list_2)
variable_df_list_all$sig = "q > 0.3"
for(i in 1:nrow(variable_df_list_all)){
  if(variable_df_list_all$q_val[i] > 0.3 ){
    variable_df_list_all$sig[i] = "q > 0.3"
  }
  if(variable_df_list_all$q_val[i] < 0.3 & variable_df_list_all$q_val[i] > 0.1 & variable_df_list_all$p_value[i] <=0.05){
    variable_df_list_all$sig[i] = "q < 0.3"
  }
  if(variable_df_list_all$q_val[i] < 0.1 ){
    variable_df_list_all$sig[i] = "q < 0.1"
  }
}


# order the factors
variable_df_list_all$Variable = factor(variable_df_list_all$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM_blast_percent', 'PB_blast_percent', 'Age'))

# plot the discrete results
p1 = ggplot(variable_df_list_all, aes(x = Mutated_Gene, y = effect_size, label = p_value)) +
  geom_hline(yintercept=0, linetype = "dashed", color = "black") +
  theme_cowplot(font_size = 10) +
  geom_pointrange(size = 0.75, stat = "identity", 
                  aes(x = Mutated_Gene, ymin = CI_lower, ymax = CI_upper, y = effect_size, color = sig)) +
  scale_color_manual(values = c( "q < 0.3" = "#1F968BFF", "q < 0.1" = "#482677FF", "q > 0.3" = "grey")) +
  ylab("Effect Size (multi-mut vs. 1 mut)")+
  xlab("")+
  theme(legend.position = "right", legend.title = element_blank(),
        axis.title.y=element_blank()) +
  coord_flip() +
  xlim(rev(levels(factor(variable_df_list_all$Mutated_Gene))))

p1 = p1 + facet_grid(. ~ Variable) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=1.5, linetype="solid"))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/multi_mut_clinical_features_discrete.pdf", dpi = 300, width = 15, height = 3, units = "in") 

# vizualize all of the survival hazard ratios

forest_plot_data_list_all$sig_color = 0

for(i in 1:nrow(forest_plot_data_list_all)){
  if(forest_plot_data_list_all$log_rank_p[i] < 0.05 & forest_plot_data_list_all$HR[i] < 1){
    forest_plot_data_list_all$sig_color[i] = 1
  }
  if(forest_plot_data_list_all$log_rank_p[i] < 0.05 & forest_plot_data_list_all$HR[i] > 1){
    forest_plot_data_list_all$sig_color[i] = 2
  }
}

forest_plot_data_list_all$sig_color = as.factor(forest_plot_data_list_all$sig_color)
forest_plot_data_list_all$sig_color = as.factor(forest_plot_data_list_all$sig_color)
forest_plot_data_list_all = subset(forest_plot_data_list_all, forest_plot_data_list_all$Upper_CI != "Inf")
forest_plot_data_list_all$categories <- reorder(forest_plot_data_list_all$categories, forest_plot_data_list_all$HR)


forest_plot_data_list_all$p_text = NA
for(i in 1:nrow(forest_plot_data_list_all)){
  if(forest_plot_data_list_all$log_rank_p[i] < 0.05){
    forest_plot_data_list_all$p_text[i] = forest_plot_data_list_all$log_rank_p[i]
  }
}
forest_plot_data_list_all$q_text = NA
for(i in 1:nrow(forest_plot_data_list_all)){
  if(forest_plot_data_list_all$q_val[i] < 0.5){
    forest_plot_data_list_all$q_text[i] = forest_plot_data_list_all$q_val[i]
  }
}

forest_plot_data_list_all$p_text = as.numeric(forest_plot_data_list_all$p_text)
forest_plot_data_list_all$q_text = as.numeric(forest_plot_data_list_all$q_text)

forest_plot_data_list_all$p_text = round(forest_plot_data_list_all$p_text, 3)
forest_plot_data_list_all$q_text = round(forest_plot_data_list_all$q_text, 2)

for(i in 1:nrow(forest_plot_data_list_all)){
  if(forest_plot_data_list_all$log_rank_p[i] < 0.01){
    forest_plot_data_list_all$p_text[i] = "p < 0.01"
  }
  if(forest_plot_data_list_all$log_rank_p[i] >= 0.01 & forest_plot_data_list_all$log_rank_p[i] <= 0.05){
    forest_plot_data_list_all$p_text[i] = paste("p =", paste(forest_plot_data_list_all$p_text[i]))
  }
}

forest_plot_data_list_all$p_q_text = paste(forest_plot_data_list_all$p_text, "; q = ", forest_plot_data_list_all$q_text, sep = "")

for(i in 1:nrow(forest_plot_data_list_all)){
  if(forest_plot_data_list_all$log_rank_p[i] > 0.05){
    forest_plot_data_list_all$p_q_text[i] = ""
  }
  if(forest_plot_data_list_all$log_rank_p[i] > 0.05){
    forest_plot_data_list_all$p_q_text[i] = ""
  }
}

forest_plot_data_list_all$HR = as.numeric(forest_plot_data_list_all$HR)
forest_plot_data_list_all$Lower_CI = as.numeric(forest_plot_data_list_all$Lower_CI)
forest_plot_data_list_all$Upper_CI = as.numeric(forest_plot_data_list_all$Upper_CI)

forest_plot_data_list_all$Mutated_Gene <- factor(forest_plot_data_list_all$Mutated_Gene, levels = forest_plot_data_list_all$Mutated_Gene[order(forest_plot_data_list_all$HR)])

p2 = ggplot(forest_plot_data_list_all, aes(x = reorder(Mutated_Gene, -HR), y = HR, label = p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(Mutated_Gene, HR+0.15), nudge_x = 0.5, nudge_y = 0, size = 2) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 19,
                  # fill = "white",
                  aes(x = Mutated_Gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("1" = "#1b7837", "2" = "#762a83"))+
  ylab("Hazard Ratio\n(multi-mut vs. 1 mut)")+
  # ylim(0,11.5)+
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/multi_mut_survival_forrest_plot.pdf", dpi = 300, width = 4, height = 3, units = "in")




### mutation category analysis ###
convergent_functional_category = final_data_matrix_2_alt %>%
  group_by(mutation_category) %>%
  dplyr::count(Sample, sort = TRUE) %>%
  subset(n>0) %>%
  na.omit() %>%
  distinct()

# add the mutational categoty data
convergent_clinical_merged_2 = convergent_clinical_merged %>%
  select(Sample, BM_blast_percent, PB_blast_percent, WBC, Hemoglobin, LDH, Platelet, Age, Censor, Time_to_OS) %>%
  unique()

convergent_clinical_merged_cat = left_join(convergent_clinical_merged_2, convergent_functional_category, by = "Sample")
convergent_clinical_merged_cat = convergent_clinical_merged_cat %>%
  select(Sample, BM_blast_percent, PB_blast_percent, WBC, Hemoglobin, LDH, Platelet, Age, Censor, Time_to_OS, mutation_category, n) %>%
  unique()

convergent_categories = final_data_matrix_2_alt %>%
  group_by(mutation_category) %>%
  dplyr::count(Sample, sort = TRUE) %>%
  distinct(mutation_category) %>%
  na.omit()

variable_df_list_2 = list()
forest_plot_data_list_2 = list()
options(scipen = 999)
k = 1
for(i in 1:nrow(convergent_categories)){
  sub_dat_final = convergent_clinical_merged_cat %>%
    subset(mutation_category == convergent_categories$mutation_category[i]) %>%
    distinct(Sample, BM_blast_percent, PB_blast_percent, WBC, Hemoglobin, LDH, Platelet, Age, Censor, Time_to_OS, n)
  
  sub_dat_final$n = ifelse(sub_dat_final$n>1, "2", "1")
  
  if(length(which(sub_dat_final$n != 1)) > 9){
    l = 1
    for(j in 2:8){
      # calculate a p-value and effect size for the difference in clinical feaures based on number of mutations
      # need to dynamically chage the column name for the variable of interest in order to calculate the effect size. will change back
      col_name = as.character(colnames(sub_dat_final[j])[j])
      names(sub_dat_final)[names(sub_dat_final) == col_name] <- "Variable"
      
      sub_dat_final[[j]] = as.numeric(sub_dat_final[[j]])
      
      # calculate effect sizes for clinical correlates
      p_val =  wilcox.test(sub_dat_final[[j]] ~ sub_dat_final$n, alternative = "two.sided")$p.value
      
      effect_size = cohens_d(data = sub_dat_final, Variable ~ n)$effsize
      effect_size_ci = cohens_d(data = sub_dat_final, Variable ~ n, ci = TRUE)
      
      variable_df <- data.frame(matrix(NA, nrow = 1, ncol = 6))
      names(variable_df) <- c("Variable", "Mutated_Category", "effect_size", "CI_lower", "CI_upper", "p_value")
      
      variable_df[1,1] <- paste(col_name)
      variable_df[1,2] <- convergent_categories$mutation_category[i]
      variable_df[1,3] <- effect_size
      variable_df[1,4] <- effect_size_ci$conf.low[1]
      variable_df[1,5] <- effect_size_ci$conf.high[1]
      variable_df[1,6] <- p_val
      
      # Add each list in the loop to a list of lists and rename column
      variable_df_list_1[[l]] = variable_df 
      l = l + 1 
      
      names(sub_dat_final)[names(sub_dat_final) == "Variable"] = col_name
      
    }
    variable_df_list_1_temp = do.call(rbind, variable_df_list_1)
    variable_df_list_1_temp = variable_df_list_1_temp %>%
      mutate(q_val = p.adjust(variable_df_list_1_temp$p_value))
    
    variable_df_list_2[[k]] = variable_df_list_1_temp
    k=k+1
    
    # calculate hazard ratios and confidence intervals for survival outcomes
    sub_dat_final$Censor = as.numeric(sub_dat_final$Censor)
    
    # run the Cox model
    model <- coxph( Surv(Time_to_OS, Censor) ~ n,
                    data = sub_dat_final )
    
    # extract the informative data from the survival model
    array_dat = summary(model)$conf.int[1:4]
    array_dat[5] = convergent_categories$mutation_category[i]
    
    # extract the log-rank p-value for the individual comparisons
    array_dat[6] = summary(model)$sctest[3]
    array_dat = array_dat[-2]
    
    forest_plot_data <- data.frame("Variable" = "Survival", "Mutated_Category" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = round(as.numeric(array_dat[2]),2), "Upper_CI" = round(as.numeric(array_dat[3]),2), "log_rank_p" = array_dat[5])
    
    forest_plot_data_list_2[[i]] = forest_plot_data
    
  }
  
}
# bind all of the clinical outcomes data and correct for FDR
forest_plot_data_list_2_all = do.call(rbind, forest_plot_data_list_2)
forest_plot_data_list_2_all$q_val = p.adjust(forest_plot_data_list_2_all$log_rank_p)

# bind all of the clinical correlates data
variable_df_list_2_all <- do.call(rbind, variable_df_list_2)
variable_df_list_all$sig = "q > 0.3"
for(i in 1:nrow(variable_df_list_2_all)){
  if(variable_df_list_2_all$q_val[i] > 0.3 ){
    variable_df_list_2_all$sig[i] = "q > 0.3"
  }
  if(variable_df_list_2_all$q_val[i] < 0.3 & variable_df_list_2_all$q_val[i] > 0.1 & variable_df_list_2_all$p_value[i] <= 0.05){
    variable_df_list_2_all$sig[i] = "q < 0.3"
  }
  if(variable_df_list_2_all$q_val[i] < 0.1){
    variable_df_list_2_all$sig[i] = "q < 0.1"
  }
}


# order the factors
variable_df_list_2_all$Variable = factor(variable_df_list_2_all$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM_blast_percent', 'PB_blast_percent', 'Age'))

# plot the discrete results
p3 = ggplot(variable_df_list_2_all, aes(x = Mutated_Category, y = effect_size, label = p_value)) +
  geom_hline(yintercept=0, linetype = "dashed", color = "black") +
  theme_cowplot(font_size = 10) +
  geom_pointrange(size = 0.75, stat = "identity", 
                  aes(x = Mutated_Category, ymin = CI_lower, ymax = CI_upper, y = effect_size, color = sig)) +
  scale_color_manual(values = c("q < 0.3" = "#1F968BFF", "q < 0.1" = "#482677FF", "q > 0.3" = "grey")) +
  ylab("Effect Size (multi-mut vs. 1 mut)")+
  xlab("")+
  theme(legend.position = "right", legend.title = element_blank(),
        axis.title.y=element_blank()) +
  coord_flip() +
  xlim(rev(levels(factor(variable_df_list_2_all$Mutated_Category))))

p3 = p3 + facet_grid(. ~ Variable) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=1.5, linetype="solid"))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/multi_category_mut_clinical_features_discrete.pdf", dpi = 300, width = 13, height = 2.5, units = "in") 

# vizualize all of the survival hazard ratios
forest_plot_data_list_2_all$sig_color = 0

for(i in 1:nrow(forest_plot_data_list_2_all)){
  if(forest_plot_data_list_2_all$log_rank_p[i] < 0.05 & forest_plot_data_list_2_all$HR[i] < 1){
    forest_plot_data_list_2_all$sig_color[i] = 1
  }
  if(forest_plot_data_list_2_all$q_val[i] < 0.05 & forest_plot_data_list_2_all$HR[i] > 1){
    forest_plot_data_list_2_all$sig_color[i] = 2
  }
  if(forest_plot_data_list_2_all$log_rank_p[i] < 0.05 & forest_plot_data_list_2_all$HR[i] > 1){
    forest_plot_data_list_2_all$sig_color[i] = 2
  }
}

forest_plot_data_list_2_all$sig_color = as.factor(forest_plot_data_list_2_all$sig_color)
forest_plot_data_list_2_all$sig_color = as.factor(forest_plot_data_list_2_all$sig_color)
forest_plot_data_list_2_all = subset(forest_plot_data_list_2_all, forest_plot_data_list_2_all$Upper_CI != "Inf")
forest_plot_data_list_2_all$categories <- reorder(forest_plot_data_list_2_all$categories, forest_plot_data_list_2_all$HR)

forest_plot_data_list_2_all$p_text = round(as.numeric(forest_plot_data_list_2_all$log_rank_p), 3)
forest_plot_data_list_2_all$q_text = as.numeric(round(forest_plot_data_list_2_all$q_val, 2))

for(i in 1:nrow(forest_plot_data_list_2_all)){
  if(forest_plot_data_list_2_all$log_rank_p[i] < 0.05){
    forest_plot_data_list_2_all$p_text[i] = forest_plot_data_list_2_all$p_text[i]
  }
  if(forest_plot_data_list_2_all$log_rank_p[i] < 0.01){
    forest_plot_data_list_2_all$p_text[i] = "p < 0.01"
  }
  if(forest_plot_data_list_2_all$q_val[i] < 0.01){
    forest_plot_data_list_2_all$q_text[i] = "q < 0.01"
  }
}

for(i in 1:nrow(forest_plot_data_list_2_all)){
  if(forest_plot_data_list_2_all$log_rank_p[i] >= 0.01 & forest_plot_data_list_2_all$log_rank_p[i] <= 0.05){
    forest_plot_data_list_2_all$p_text[i] = paste("p =", paste(forest_plot_data_list_2_all$p_text[i]))
  }
  if(forest_plot_data_list_2_all$q_val[i] >= 0.01 & forest_plot_data_list_2_all$q_val[i] <= 0.05){
    forest_plot_data_list_2_all$q_text[i] = paste("q =", paste(forest_plot_data_list_2_all$q_text[i]))
  }
}

forest_plot_data_list_2_all$p_q_text = paste(forest_plot_data_list_2_all$p_text, "; ", forest_plot_data_list_2_all$q_text, sep = "")

for(i in 1:nrow(forest_plot_data_list_2_all)){
  if(forest_plot_data_list_2_all$log_rank_p[i] > 0.05){
    forest_plot_data_list_2_all$p_q_text[i] = ""
  }
  if(forest_plot_data_list_2_all$log_rank_p[i] > 0.05){
    forest_plot_data_list_2_all$p_q_text[i] = ""
  }
}

forest_plot_data_list_2_all$HR = as.numeric(forest_plot_data_list_2_all$HR)
forest_plot_data_list_2_all$Lower_CI = as.numeric(forest_plot_data_list_2_all$Lower_CI)
forest_plot_data_list_2_all$Upper_CI = as.numeric(forest_plot_data_list_2_all$Upper_CI)

forest_plot_data_list_2_all$Mutated_Category <- factor(forest_plot_data_list_2_all$Mutated_Category, levels = forest_plot_data_list_2_all$Mutated_Category[order(forest_plot_data_list_2_all$HR)])

p4 = ggplot(forest_plot_data_list_2_all, aes(x = reorder(Mutated_Category, -HR), y = HR, label = p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(Mutated_Category, HR+0.15), nudge_x = 0.5, nudge_y = 0, size = 2) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 19,
                  aes(x = Mutated_Category, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("1" = "#1b7837", "2" = "#762a83"))+
  ylab("Hazard Ratio\n(multi-mut vs. 1 mut)")+
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/multi_category_mut_survival_forrest_plot.pdf", dpi = 300, width = 4, height = 3, units = "in")

# arrange all of the plots together
library(patchwork)
p1 + p2 + p3 + p4 + plot_layout(widths = c(5, 1.5), heights = c(1,.75))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_2/Supplimental/convergent_panel_2.pdf", dpi = 300, width = 15, height = 5.75, units = "in")
