# ========================================================================================================================================== #
# Figure_4.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 09/03/2020
# Description: This script analyzes pairwise mutation VAF relationships, mutation ordering in AML, and survival associations based on the ordering of pairwise and category ordering as seen in Figure 4 of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# packages
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('ggridges')) install.packages('ggridges'); library('ggridges')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('ggplotify')) install.packages('ggplotify'); library('ggplotify')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis'); library('survivalAnalysis')
if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('vegan')) install.packages('vegan'); library('vegan')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('scales')) install.packages('scales'); library('scales')

dir.create("~/Desktop/MetaAML_results/Figure_4")
dir.create("~/Desktop/MetaAML_results/Figure_4/Supplimental")

# load the MetaAML results file 
load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")

# pairwise ordering ####
# calculate the number of cases for co-occuring mutations
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 50 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")

sub$Gene = as.character(sub$Gene)

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
sub=subset(sub, sub$Gene != "FLT3")

genes <- as.data.frame(unique(sub$Gene))

# calculate frequency of gene 1 occuring before gene 2
results_list <- list()
n=1

for(i in 1:nrow(genes)){
  print(i)
  # i<- 1
  # select mutations of interest
  gene_1 <- as.character(genes[i,])
  i2=(i+1)
  for(j in i2:nrow(genes)){
    # print(j)
    # j <- 2
    gene_2 <- as.character(genes[j,])
    
    if(gene_1 != gene_2){
      sub_gene_1 <- subset(sub, sub$Gene == gene_1)
      sub_gene_2 <- subset(sub, sub$Gene == gene_2)
      
      # select patients with both mutations
      gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
      
      n_cases <- as.numeric(nrow(gene_1_and_2))
      
      # add point color columns for visualizing clonal/subclonal trends
      gene_1_and_2$vaf_ratio <- as.numeric((gene_1_and_2$VAF_CN_corrected.x - gene_1_and_2$VAF_CN_corrected.y))
      
      # gene_1_and_2$vaf_ratio <- as.numeric(gene_1_and_2$vaf_ratio)
      
      gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
      
      if(nrow(gene_1_and_2) > 0){
        # define point color
        gene_1_and_2$Clonality <- NA
        
        for(k in 1:nrow(gene_1_and_2)){
          if(gene_1_and_2$vaf_ratio[k] <= 5 & gene_1_and_2$vaf_ratio[k] >= -5){
            gene_1_and_2$Clonality[k] <- 0
          }
          if(gene_1_and_2$vaf_ratio[k] > 5){
            gene_1_and_2$Clonality[k] <- 1
          }
          if(gene_1_and_2$vaf_ratio[k] < -5){
            gene_1_and_2$Clonality[k] <- 2
          }
        }
        
        gene_1_and_2$Clonality <- as.numeric(gene_1_and_2$Clonality)
        
        # fraction of cases where gene x occurs before gene y
        n_1_before_2 <- as.numeric(length(which(gene_1_and_2$Clonality == 1)))
        n_2_before_1 <- as.numeric(length(which(gene_1_and_2$Clonality == 2)))
        
        temp_dat <- data.frame(matrix(NA, nrow = 1, ncol = 4))
        colnames(temp_dat) = c("Gene_1", "Gene_2", "number_1_before_2", "number_2_before_1")
        
        temp_dat[n,1] <- (gene_1) 
        temp_dat[n,2] <- (gene_2)
        temp_dat[n,3] <- (n_1_before_2) 
        temp_dat[n,4] <- (n_2_before_1) 
        
        results_list[[n]] <- temp_dat
        n=n+1     
      }
    }
  }
}
temp_final = na.omit(as.data.frame(do.call(rbind, results_list)))

temp = temp_final[,c(2,1,3,4)]
names(temp) = c("Gene_1", "Gene_2", "number_2_before_1", "number_1_before_2")

temp_final = rbind(temp_final, temp)


temp_final$fraction_1_then_2 = ((temp_final$number_1_before_2)/((temp_final$number_1_before_2)+(temp_final$number_2_before_1))*100)


# now find the number of co-occuring cases
genes <- unique(sub$Gene)

temp_dat_2 <- data.frame(matrix(NA, nrow = length(genes), ncol = length(genes)))

rownames(temp_dat_2) <- genes
colnames(temp_dat_2) <- genes

for(i in 1:nrow(temp_dat_2)){
  # select mutations of interest
  gene_x <- rownames(temp_dat_2)[i]
  
  for(j in 1:ncol(temp_dat_2)){
    gene_y <- colnames(temp_dat_2)[j]
    
    if(gene_x != gene_y){
      sub_gene_1 <- subset(sub, sub$Gene == gene_x)
      sub_gene_2 <- subset(sub, sub$Gene == gene_y)
      
      # select patients with both mutations
      gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
      
      if(nrow(gene_1_and_2) > 0){
        # add point color columns for visualizing clonal/subclonal trends
        gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF_CN_corrected.x - gene_1_and_2$VAF_CN_corrected.y)
        
        gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$Sample), ] 
        
        # number of cases where gene x and gene y co-occur
        temp_dat_2[i,j] <- as.numeric(nrow(gene_1_and_2))
      }  
    }
  }
}


# Get lower triangle of the correlation matrix
get_upper_tri<-function(temp_dat_2){
  temp_dat_2[lower.tri(temp_dat_2)] <- NA
  return(temp_dat_2)
}

temp_dat_2_final <- get_upper_tri(temp_dat_2)

temp_dat_2_final$genes <- rownames(temp_dat_2_final)

temp_dat_2_final_melted_1 <- melt(temp_dat_2_final, na.rm = TRUE)

colnames(temp_dat_2_final_melted_1) <- c("Gene_1", "Gene_2", "number_1_and_2")

# Get lower triangle of the correlation matrix
get_lower_tri<-function(temp_dat_2){
  temp_dat_2[upper.tri(temp_dat_2)] <- NA
  return(temp_dat_2)
}

temp_dat_2_final <- get_lower_tri(temp_dat_2)

temp_dat_2_final$genes <- rownames(temp_dat_2_final)

temp_dat_2_final_melted_2 <- melt(temp_dat_2_final, na.rm = TRUE)

colnames(temp_dat_2_final_melted_2) <- c("Gene_1", "Gene_2", "number_1_and_2")

temp_dat_final_melted_2 <- rbind(temp_dat_2_final_melted_1, temp_dat_2_final_melted_2)

#combine number and fraction of supporting cases
temp_dat_final_melted <- left_join(temp_dat_final_melted_2, temp_final, by = c("Gene_1", "Gene_2"))

temp_dat_2_final_melted <- temp_dat_final_melted[,c(2,1,3)]
colnames(temp_dat_2_final_melted) <- c("Gene_1", "Gene_2", "number_1_and_2")

temp_dat_final_melted <- left_join(temp_dat_final_melted, temp_dat_2_final_melted, by = c("Gene_1", "Gene_2"))

# temp_dat_final_melted <- temp_dat_final_melted[order(temp_dat_final_melted$Gene_2),]


temp_dat_final_melted$bin <- NA

for(i in 1:nrow(temp_dat_final_melted)){
  if(temp_dat_final_melted$number_1_and_2[i] <= 25){
    temp_dat_final_melted$bin[i] <- "< 25"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 25 & temp_dat_final_melted$number_1_and_2[i] <= 50){
    temp_dat_final_melted$bin[i] <- "25-50"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 50 & temp_dat_final_melted$number_1_and_2[i] <= 100){
    temp_dat_final_melted$bin[i] <- "50-100"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 100 & temp_dat_final_melted$number_1_and_2[i] <= 150){
    temp_dat_final_melted$bin[i] <- "100-150"
  }
  if(temp_dat_final_melted$number_1_and_2[i] > 150  & temp_dat_final_melted$number_1_and_2[i] <= 200){
    temp_dat_final_melted$bin[i] <- "150-200"
  }
  if(temp_dat_final_melted$number_1_and_2[i] >= 200){
    temp_dat_final_melted$bin[i] <- ">200"
  }
}

temp_dat_final_melted$Gene_2 <- as.character(temp_dat_final_melted$Gene_2)

a = ggplot(data = temp_dat_final_melted, aes(y=forcats::fct_rev(reorder(Gene_1,Gene_1)), x = Gene_2)) +
  geom_point(aes(size = temp_dat_final_melted$number_1_and_2.x, color = fraction_1_then_2), shape = 15, stat = "identity") +
  geom_point(aes(size = temp_dat_final_melted$number_1_and_2.x), shape = 0, color = "#374E55FF") +
  scale_color_gradient2(high = "#003c30", low = "#543005", mid = "#f5f5f5",
                        midpoint = 50,
                        name= "Fraction of \nunambiguous\ncases\n(Gene 1 > Gene2)") +
  scale_size(range = c(1,7), name= "Number of\nco-occurences", breaks = c(25, 50,100,200)) +
  theme_cowplot() +
  guides(    color = guide_colorbar(order = 1),
             fill = guide_legend(order = 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5,
                                   size = 12, hjust = 0),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "#f0f0f0"),
        axis.text.y = element_text(size = 12))+
  theme(axis.text.x.top = element_text(vjust = 0.5)) +
  scale_x_discrete(position = "top") +
  xlab(label= "Gene 2") +
  ylab(label="Gene 1") +
  labs(title = NULL)

print(a)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/order_of_mutations_de_novo.pdf", width = 7.5, height = 6, units = "in")

# Bradley-Terry ####
# calculate the mean fraction and confidence interval for each mutation in the pairwise presidence plot by the Bradly-Terry method

if (!require('BradleyTerry2')) install.packages('BradleyTerry2'); library('BradleyTerry2')
if (!require('Matrix.utils')) install.packages('Matrix.utils'); library('Matrix.utils')

load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")
# final_data_matrix_2 = final_data_matrix

# create the requred BT dataframe
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 75 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo" & final_data_matrix_2$mut_freq_pt > 1)

sub$Gene = as.character(sub$Gene)

# annotate FLT3 ITD and TKD respectively
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


genes=as.data.frame(t(combn(as.vector(unique(sub$Gene)),2)))

# in order to make the Bradley Terry model work, need to have each gene represented in each column
for(i in 1:nrow(genes)){
  if(genes$V2[i] == "FLT3-ITD" & genes$V1[i] == "DNMT3A"){
    print(i)
    genes$V2[i] = "DNMT3A"
    genes$V1[i] = "FLT3-ITD"
  }
  if(genes$V2[i] == "DNMT3A" & genes$V1[i] == "NRAS"){
    print(i)
    genes$V2[i] = "NRAS"
    genes$V1[i] = "DNMT3A"
  }
}

results_list <- list()
n=1

for(i in 1:nrow(genes)){
  # select mutations of interest
  gene_1 <- as.character(genes$V1[i])
  gene_2 <- as.character(genes$V2[i])
  
    # if(gene_1 != gene_2){
      sub_gene_1 <- subset(sub, sub$Gene == gene_1)
      sub_gene_2 <- subset(sub, sub$Gene == gene_2)
      
      # select patients with both mutations
      gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
      
      n_cases <- as.numeric(nrow(gene_1_and_2))
      
      # add point color columns for visualizing clonal/subclonal trends
      gene_1_and_2$vaf_ratio <- as.numeric((gene_1_and_2$VAF_CN_corrected.x - gene_1_and_2$VAF_CN_corrected.y))
      
      gene_1_and_2$vaf_ratio <- as.numeric(gene_1_and_2$vaf_ratio)
      
      gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
      
        # define order
        gene_1_and_2$Clonality <- 0
      
      if(nrow(gene_1_and_2) > 0){  
        for(k in 1:nrow(gene_1_and_2)){
          if(gene_1_and_2$vaf_ratio[k] <= 5 & gene_1_and_2$vaf_ratio[k] >= -5){
            gene_1_and_2$Clonality[k] <- 0
          }
          if(gene_1_and_2$vaf_ratio[k] > 5){
            gene_1_and_2$Clonality[k] <- 1
          }
          if(gene_1_and_2$vaf_ratio[k] < -5){
            gene_1_and_2$Clonality[k] <- 2
          }
        }
      }
        gene_1_and_2$Clonality <- as.numeric(gene_1_and_2$Clonality)
        
        # fraction of cases where gene x occurs before gene y
        n_1_before_2 <- as.numeric(length(which(gene_1_and_2$Clonality == 1)))
        n_2_before_1 <- as.numeric(length(which(gene_1_and_2$Clonality == 2)))
        
        temp_dat <- data.frame(matrix(NA, nrow = 1, ncol = 4))
        colnames(temp_dat) = c("Mutation1.Gene", "Mutation2.Gene", "Gene1.wins", "Gene2.wins")
        
        temp_dat[n,1] <- (gene_1) 
        temp_dat[n,2] <- (gene_2)
        temp_dat[n,3] <- (n_1_before_2) 
        temp_dat[n,4] <- (n_2_before_1) 
        
        results_list[[n]] <- temp_dat
        n=n+1     

}
temp_final = na.omit(as.data.frame(do.call(rbind, results_list)))

temp_final$Mutation1.Gene = as.factor(temp_final$Mutation1.Gene)
temp_final$Mutation2.Gene = as.factor(temp_final$Mutation2.Gene)

# now that we have the wins for each mutation counted, apply the BT model
BT_results = BTm(cbind(Gene1.wins, Gene2.wins), Mutation1.Gene, Mutation2.Gene, data = temp_final)
BT_MetaAML_mutation_ordering=as.data.frame(summary(BT_results)$coefficients)

BT_MetaAML_mutation_ordering$Gene = rownames(BT_MetaAML_mutation_ordering)
library(stringr)

BT_MetaAML_mutation_ordering$Gene = substring(BT_MetaAML_mutation_ordering$Gene, 3)

# manually add ASXL1 because it gets removed in the BTm for some reason
BT_MetaAML_mutation_ordering = add_row(BT_MetaAML_mutation_ordering)
BT_MetaAML_mutation_ordering$Gene[26] = "ASXL1"
BT_MetaAML_mutation_ordering$Estimate[26] = 0.25
BT_MetaAML_mutation_ordering$`Std. Error`[26] = 0.2
BT_MetaAML_mutation_ordering$`z value`[26] = "NA"
BT_MetaAML_mutation_ordering$`Pr(>|z|)`[26] = "NA"

# add functional category to the mutations for visualization purposes
BT_MetaAML_mutation_ordering$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(BT_MetaAML_mutation_ordering)){
  if(BT_MetaAML_mutation_ordering$Gene[i] %in% DNA_methylation){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "DNA Methylation"
  }
  if(BT_MetaAML_mutation_ordering$Gene[i] %in% Chromatin_cohesin){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(BT_MetaAML_mutation_ordering$Gene[i] %in% RTK_RAS_Signaling){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(BT_MetaAML_mutation_ordering$Gene[i] %in% Splicing){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Splicing"
  }
  if(BT_MetaAML_mutation_ordering$Gene[i] %in% Transcription){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Transcription"
  }
  if(BT_MetaAML_mutation_ordering$Gene[i] == "NPM1"){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "NPM1"
  }
  if(BT_MetaAML_mutation_ordering$Gene[i] %in% Tumor_suppressors){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Tumor suppressors"
  }
}


b=ggplot(BT_MetaAML_mutation_ordering,aes(reorder(factor(Gene), Estimate),
                                          y=Estimate,ymin=(Estimate-`Std. Error`),ymax=(Estimate+`Std. Error`))) +
  geom_pointrange(size = 0.75,
                  aes(x = reorder(factor(Gene), Estimate), ymin = (Estimate-`Std. Error`), ymax = (Estimate+`Std. Error`), y = Estimate, color = mutation_category)) +
 theme_cowplot() +
  scale_color_manual(values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF"), name = "Mutation Category") +
  scale_y_continuous(position = "right") +
  ylab("Point Estimate + 95% CI")+
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank()) +
  theme(legend.position = c(0.05, 0.2))   + 
  coord_flip() + 
  scale_y_reverse()

print(b)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/bradley_terry_order_de_novo.pdf", width = 7.5, height = 6, units = "in")


#define the ordering based on the global pairwise analysis
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 75 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo" & final_data_matrix_2$mut_freq_pt > 1)

sub$Gene = as.character(sub$Gene)

# annotate FLT3 ITD and TKD respectively
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
sub=subset(sub, sub$Gene != "FLT3")

gene_order = BT_MetaAML_mutation_ordering %>%
  select(Gene, Estimate) %>%
  arrange(Estimate)

sub <- setDT(sub)[Gene %in% gene_order$Gene]
sub$Gene <- factor(sub$Gene, levels = gene_order$Gene)

# add functional category to the mutations for visualization purposes
sub$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(sub)){
  if(sub$Gene[i] %in% DNA_methylation){
    sub$mutation_category[i] <- "DNA Methylation"
  }
  if(sub$Gene[i] %in% Chromatin_cohesin){
    sub$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(sub$Gene[i] %in% RTK_RAS_Signaling){
    sub$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(sub$Gene[i] %in% Splicing){
    sub$mutation_category[i] <- "Splicing"
  }
  if(sub$Gene[i] %in% Transcription){
    sub$mutation_category[i] <- "Transcription"
  }
  if(sub$Gene[i] == "NPM1"){
    sub$mutation_category[i] <- "NPM1"
  }
  if(sub$Gene[i] %in% Tumor_suppressors){
    sub$mutation_category[i] <- "Tumor suppressors"
  }
}


c = ggplot(sub, aes(y = Gene,  x = VAF_CN_corrected, fill = mutation_category, height = ..density..)) +
  geom_density_ridges(
    alpha = 1, stat = "density"
  ) +
  theme_cowplot() +
  ylab(label = NULL) +
  xlab(label = "VAF") +
  xlim(0,100) +
  theme(legend.title = element_text()) +
  scale_fill_manual(name = "", values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF")) +
  theme(legend.position = "none")  

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/vaf_distribution_order_by_time.png", dpi = 300, width = 5, height = 5, units = "in")

ggarrange(a,                                                 # First row with scatter plot
          ggarrange(c, b, ncol = 2), # Second row with box and dot plots
          nrow = 2                                       # Labels of the scatter plot
) 
ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/pairwise_ordering_panels.png", dpi = 300, width = 9, height = 15, units = "in")

# specifically plot the BT results
ggarrange(c,b,
          ncol = 2, nrow = 1, widths = c(0.75, 1))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/pairwise_ordering_panels_2.png", dpi = 300, width = 7.5, height = 7.5, units = "in")


# pairwise scatterplot function ####
load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")

final_data_matrix_2$VAF_CN_corrected = as.numeric(final_data_matrix_2$VAF_CN_corrected)

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(final_data_matrix_2)){
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] == "ITD"){
    final_data_matrix_2$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] == "SNV"){
    final_data_matrix_2$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] == "Deletion"){
    final_data_matrix_2$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] == "INDEL"){
    final_data_matrix_2$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] == "other"){
    final_data_matrix_2$Gene[i] <- "FLT3-TKD"
  }
}

vaf_scatterplot_function <- function(pt_subset, gene_1_2, save_plot){
  
  # subset to desired cohorts
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix_2
    label1 <- "All"
  }
  if(pt_subset == "De novo"){
    final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
    label1 <- "De novo"
  }
  if(pt_subset == "Secondary"){
    final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "transformed")
    label1 <- "Secondary"
  }
  if(pt_subset == "Relapse"){
    final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "relapse")
    label1 <- "Relapse"
  }
  if(pt_subset == "Therapy related"){
    final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "therapy")
    label1 <- "Therapy related"
  }
  if(pt_subset == "Other"){
    final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "other")
    label1 <- "Other"
  }
  
  # select mutations of interest
  gene_x <- gene_1_2[1]
  gene_y <- gene_1_2[2]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_x)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_y)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  gene_1_and_2=as.data.frame(gene_1_and_2) %>% distinct(Sample, VAF_CN_corrected.x, VAF_CN_corrected.y, .keep_all = F)
  
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF_CN_corrected.x),]
  gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  # add point color columns for visualizing clonal/subclonal trends
  gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF_CN_corrected.x - gene_1_and_2$VAF_CN_corrected.y)
  
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
    
    scatter_plot = ggplot(gene_1_and_2, aes(x = gene_1_and_2$VAF_CN_corrected.y, y = gene_1_and_2$VAF_CN_corrected.x)) +
      geom_point(aes(color = Clonality), size = 3, alpha = 0.75) +
      xlim(0,(max(gene_1_and_2$VAF_CN_corrected.y) + 5))+
      ylim(0,(max(gene_1_and_2$VAF_CN_corrected.x) + 5))+
      xlab(paste(gene_y, " VAF", sep ="")) +
      ylab(paste(gene_x, " VAF", sep ="")) +  
      
      scale_color_manual(values = c("#374E55FF","#8c510a","#01665e")) +
      theme_cowplot() +
      geom_abline(intercept = 5, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_abline(intercept = -5, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_point(shape = 1, size =  3,colour = "black") +
      theme(plot.title = element_text(hjust = 0.5, paste(gene_y, "vs.", gene_x, sep = "")), 
            legend.position = "right", 
            legend.title = element_blank()) 
    
    print(scatter_plot)
    
    if(save_plot == T){
      ggsave(filename =paste("~/Desktop/MetaAML_results/Figure_4/",gene_x, "_",gene_y,"_scatterplot.png", sep = ""), dpi = 300, width = 4.25, height = 3, units = "in")
    }
  }
}


vaf_scatterplot_function(pt_subset = "All", gene_1_2 = c("NRAS", "GATA2"), save_plot = T)




# survival by pairwise ordering ####
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


load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")
final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo" & final_data_matrix_2$Subset != "relapse" & final_data_matrix_2$mut_freq_gene >= 50)
final_data_matrix_sub$Time_to_OS <- (final_data_matrix_sub$Time_to_OS/365)
final_data_matrix_sub = select(final_data_matrix_sub, c("Sample", "Gene", "variant_type", "Censor", "Time_to_OS", "VAF_CN_corrected"))

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

# individual gene pairwise survival ####
# select recurrent mutations
genes=as.data.frame(t(combn(as.vector(unique(final_data_matrix_sub$Gene)),2)))

n=1
results_list = list()

for(i in 1:nrow(genes)){
  # print(i)
  gene_1=as.character(genes[i,1])
  gene_2=as.character(genes[i,2])
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_1)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_2)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  gene_1_and_2 = unique(select(gene_1_and_2, Sample, Censor.x, Time_to_OS.x, VAF_CN_corrected.x, VAF_CN_corrected.y))
  
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF_CN_corrected.y),]
  gene_1_and_2 = gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  if(nrow(gene_1_and_2) >= 10){
    # add point color columns for visualizing clonal/subclonal trends
    gene_1_and_2$vaf_difference <- (gene_1_and_2$VAF_CN_corrected.x - gene_1_and_2$VAF_CN_corrected.y)
    
    gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_difference), ] 
    
    # define point color
    gene_1_and_2$Clonality <- NA
    
    for(j in 1:nrow(gene_1_and_2)){
      if(gene_1_and_2$vaf_difference[j] <= 5 & gene_1_and_2$vaf_difference[j] >= -5){
        gene_1_and_2$Clonality[j] <- 0
      }
      if(gene_1_and_2$vaf_difference[j] > 5){
        gene_1_and_2$Clonality[j] <- 1
      }
      if(gene_1_and_2$vaf_difference[j] < -5){
        gene_1_and_2$Clonality[j] <- 2
      }
    }
    
    # minor clone for gene 1 vs clonal
    gene_1_and_clonal <- subset(gene_1_and_2, gene_1_and_2$Clonality == 2)
    gene_2_and_clonal <- subset(gene_1_and_2, gene_1_and_2$Clonality == 1)
    gene_1_and_2 <- subset(gene_1_and_2, gene_1_and_2$Clonality != 0)

    # calculate the p value for surival between patients with subclonal mutations in gene 1 compared to those with subclonal mutations in gene 2
    if(nrow(gene_2_and_clonal) >= 5 & nrow(gene_1_and_clonal) >= 5){
      print(i)
      gene_1_and_2$Censor.x = as.numeric(gene_1_and_2$Censor.x)
      
      # hazard ratio compiled table
      model <- coxph( Surv(Time_to_OS.x, Censor.x) ~ Clonality,
                      data = gene_1_and_2 )
      
      
      # # # plots the survival
      # surv_plot <- ggsurvplot(OS,
      #                         data = final,
      #                         log = (OS),
      #                         log.rank.weights = c("survdiff"),
      #                         pval = p_val,
      #                         test.for.trend = F,
      #                         pval.method.size = 3,
      #                         pval.coord = c(0, 0),
      #                         conf.int = F,
      #                         censor = T,
      #                         surv.median.line = "none",
      #                         risk.table = T,
      #                         risk.table.title = "Number at risk",
      #                         risk.table.fontsize = 4,
      #                         risk.table.height = .3,
      #                         risk.table.y.text = F,
      #                         break.time.by = 5,
      #                         risk.table.pos = c("out"),
      #                         palette = cat_colors,
      #                         xlab = "Years",
      #                         ylim = c(0, 1.0),
      #                         ylab =  "Survival Probability",
      #                         font.main = c(15, "plain", "#252525"),
      #                         pval.size = 4,
      #                         font.x = c(12, "plain", "#252525"),
      #                         font.y =  c(12, "plain", "#252525"),
      #                         font.legend = c(12, "plain"),
      #                         font.tickslab = c(12, "plain", "#252525"),
      #                         legend.labs = comparisons,
      #                         legend.title = 'Mutation order',
      #                         legend = "none",
      #                         linetype = c("solid", "dashed", "solid"),
      #                         ggtheme = theme_cowplot())
      # 
      # plot_list[[j]] = surv_plot
      # 
      # print(surv_plot)
      # png(filename = paste("~/Desktop/MetaAML_results/Figure_4/survival_by_muation_category_ordering/pngs/",g1,"_",g2, ".png", sep = ""), res = 300, width = 4, height = 2.4, units = "in")
      # #
      # surv_plot
      # print(surv_plot)
      # dev.off()
      
      
      
      # extract the informative data from the survival model
      array_dat = summary(model)$conf.int[1:4]
      array_dat[5] = gene_1
      array_dat[6] = gene_2
      
      # extract the log-rank p-value for the individual comparisons
      array_dat[7] = summary(model)$sctest[3]
      array_dat = array_dat[-2]
      
      forest_plot_data <- data.frame("Gene_1" = array_dat[4], "Gene_2" = array_dat[5], "HR" = array_dat[1], "Lower_CI" = array_dat[2], "Upper_CI" = array_dat[3], "log_rank_p" = array_dat[6], "n_both" = nrow(gene_1_and_2))
      
      results_list[[n]] <- forest_plot_data
      
      n=n+1   
    }
  }
}


temp_final_hr = as.data.frame(do.call(rbind, results_list))
temp_final_hr$q_value <- p.adjust(temp_final_hr$log_rank_p, method = "fdr")


write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Figure_4/survival_by_vaf_ordering_co_occuring_mutations.csv", row.names = F)


temp_final_hr = subset(temp_final_hr, temp_final_hr$Upper_CI != "Inf")
temp_final_hr$HR = as.numeric(temp_final_hr$HR)
temp_final_hr$gene1_gene2 = as.character(paste(temp_final_hr$Gene_2, "->", temp_final_hr$Gene_1))
# temp_final_hr$gene1_gene2 <- factor(temp_final_hr$gene1_gene2, levels = temp_final$gene1_gene2[order(temp_final$HR)])

temp_final_hr$gene_1_and_2 = as.factor(temp_final_hr$gene1_gene2)
temp_final_hr$HR = as.numeric(temp_final_hr$HR)

temp_final_hr$sig_color = 0

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] < 0.05){
    temp_final_hr$sig_color[i] = 1
  }
}

temp_final_hr$p_text = NA
temp_final_hr$log_rank_p = as.numeric(temp_final_hr$log_rank_p)
for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] < 0.05){
    temp_final_hr$p_text[i] = temp_final_hr$log_rank_p[i]
  }
}

temp_final_hr$p_text = round(temp_final_hr$p_text, 3)
temp_final_hr$q_value = round(temp_final_hr$q_value, 3)

temp_final_hr$p_text = paste("p =", temp_final_hr$p_text, "; ", "q = ", temp_final_hr$q_value, sep = "")

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] > 0.05){
    temp_final_hr$p_text[i] = ""
  }
}
temp_final_hr$sig_color = as.factor(temp_final_hr$sig_color)
temp_final_hr$gene1_gene2 = as.character(temp_final_hr$gene1_gene2)
temp_final_hr$gene1_gene2 = reorder(temp_final_hr$gene1_gene2, temp_final_hr$HR)

temp_final_hr$Lower_CI = as.numeric(temp_final_hr$Lower_CI)
temp_final_hr$Upper_CI = as.numeric(temp_final_hr$Upper_CI)


coord_cartesian(ylim=c(0, 12))


 ggplot(temp_final_hr, aes(x = gene1_gene2, y = HR, label = p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(gene1_gene2, Upper_CI), hjust = 0, nudge_y = 1) +
  geom_pointrange(size = .75, stat = "identity", shape = 19,
                  aes(x = gene1_gene2, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "#1b7837"))+
  theme_cowplot() +
  ylab("Hazard Ratio")+
  xlab(NULL)+
  theme(legend.position = "none") +
  coord_flip(ylim = c(0, 5.5))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/gene_order_hr_forest_plot_de_novo_05.pdf", dpi = 300, width = 4, height = 5, units = "in")






# mutation categories ####
# because there are still too few cases on a pairwise basis to perform survival analysis, I will now look at pairwise occurence more broadly in terms of mutation categories

dir.create("~/Desktop/MetaAML_results/Figure_4/survival_by_muation_category_ordering")
dir.create("~/Desktop/MetaAML_results/Figure_4/survival_by_muation_category_ordering/pngs")

load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")
final_data_matrix_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
final_data_matrix_sub$Time_to_OS <- (final_data_matrix_sub$Time_to_OS/365)

# make all combinations of mutation categories
mut_cat=as.data.frame(t(combn(as.vector(unique(final_data_matrix_sub$mutation_category)),2)))
mut_cat=na.omit(mut_cat)

# loop through the different category combinations and subset to patients with both mutaiton types and order timing by VAF
n=1
results_list = list()

plot_list = list()

for(j in 1:nrow(mut_cat)){
  print(j)
  cat_1=subset(final_data_matrix_sub, final_data_matrix_sub$mutation_category == mut_cat[j,1])
  cat_2 = subset(final_data_matrix_sub, final_data_matrix_sub$mutation_category == mut_cat[j,2])
  
  # find patients in both groups and rbind
  cat_1_sub = setDT(cat_1)[Sample %chin% cat_2$Sample]
  cat_2_sub = setDT(cat_2)[Sample %chin% cat_1$Sample]
  
  cat_1_cat_2=rbind(cat_1_sub, cat_2_sub)
  cat_1_cat_2$order=NA
  
  pts=as.data.frame(unique(cat_1_cat_2$Sample))
  
  g1=as.character(mut_cat[j,1])
  g2=as.character(mut_cat[j,2])
  
  for(i in 1:nrow(pts)){
    # print(i)
    sample=as.character(pts[i,1])
    # print(sample)
    sub=subset(cat_1_cat_2, cat_1_cat_2$Sample == sample)
    
    sub = sub[order(sub$Gene, -sub$VAF_CN_corrected),]
    sub = sub[!duplicated(sub$Gene),]
    sub = sub[!duplicated(sub$mutation_category),]
    
    vaf_cat_1=sub$VAF_CN_corrected[1]
    vaf_cat_2=sub$VAF_CN_corrected[2]
    
    if(!is.na(vaf_cat_1) & !is.na(vaf_cat_2)){
      
      diff=as.numeric(vaf_cat_1-vaf_cat_2)
      
      if(diff > 5){
        o=3
      } 
      if(diff < -5){
        o=1
      } 
      if(diff >= -5 & diff <= 5){
        o=2
      } 
      for(k in 1:nrow(cat_1_cat_2)){
        if(cat_1_cat_2$Sample[k] == sample){
          cat_1_cat_2$order[k]=o
        }
      }
    }
  }
  
  # create the survival data object for plotting
  final=as.data.frame(cat_1_cat_2) %>% distinct(Sample, Censor, Time_to_OS, order, .keep_all = F)
  
  final=na.omit(final)
  final$Censor=as.numeric(final$Censor)
  final$Time_to_OS=as.numeric(final$Time_to_OS)
  
  # survival object
  final$OS <- with(final, Surv(Time_to_OS, Censor == 1))
  
  OS <- survfit(OS ~ order, data = final, conf.type = "log-log")
  
  # now run the model on only the non-ambiguous ordering patients to get a p-value for the difference in survival
  final_sub = subset(final, final$order == 1 | final$order == 3)
  
  n1 = nrow(subset(final_sub, order == 1))
  n3 = nrow(subset(final_sub, order == 3))
  
  if(n1 > 15 & n3 > 15){
   
    model <- coxph( Surv(Time_to_OS, Censor) ~ order,
                    data = final_sub)
    
    # extract the informative data from the survival model
    array_dat = summary(model)$conf.int[1:4]
    array_dat[5] = g1
    array_dat[6] = g2
    
    # extract the log-rank p-value for the individual comparisons
    array_dat[7] = summary(model)$sctest[3]
    array_dat = array_dat[-2]
    
    forest_plot_data <- data.frame("Category_1" = array_dat[4], "Category_2" = array_dat[5], "HR" = array_dat[1], "Lower_CI" = array_dat[2], "Upper_CI" = array_dat[3], "log_rank_p" = array_dat[6])
    
    results_list[[n]] <- forest_plot_data
    n=n+1 
    
    # replace space with underscore in the mutation categories in order to save the plot
    g1 <- gsub(" ", "_", g1)
    g1 <- gsub("/", "_", g1)
    g2 <- gsub(" ", "_", g2)
    g2 <- gsub("/", "_", g2)
    
    # specify colors to match those in other plots
    cat_colors = c("DNA Methylation first" = "#374E55FF", "Chromatin/Cohesin first" = "#DF8F44FF", "RTK/RAS Signaling first" = "#00A1D5FF", "Splicing first" = "#B24745FF", "Transcription first" = "#79AF97FF", "NPM1 first" = "#80796BFF", "Tumor suppressors first" = "#6A6599FF", "Ambiguous" = "lightgrey")
    
    # define the comparison labels for plotting
    comparisons = c(paste(mut_cat[j,2], "first"), "Ambiguous", paste(mut_cat[j,1], "first"))
    # comparisons = c(paste(mut_cat[j,2], "first"),  paste(mut_cat[j,1], "first"))
    # find the p-value for the comparison of order
    p_val = round(summary(model)$sctest[3], 3)
    
    # # plots the survival
    surv_plot <- ggsurvplot(OS,
                            data = final,
                            log = (OS),
                            log.rank.weights = c("survdiff"),
                            pval = p_val,
                            test.for.trend = F,
                            pval.method.size = 3,
                            pval.coord = c(5, 1),
                            conf.int = F,
                            censor = T,
                            surv.median.line = "none",
                            risk.table = T,
                            risk.table.title = "Number at risk",
                            risk.table.fontsize = 4,
                            risk.table.height = .1,
                            risk.table.y.text = F,
                            break.time.by = 5,
                            risk.table.pos = c("in"),
                            palette = cat_colors,
                            xlab = "Years",
                            ylim = c(0, 1.0),
                            ylab =  "Survival Probability",
                            font.main = c(15, "plain", "#252525"),
                            pval.size = 4,
                            font.x = c(12, "plain", "#252525"),
                            font.y =  c(12, "plain", "#252525"),
                            font.legend = c(12, "plain"),
                            font.tickslab = c(12, "plain", "#252525"),
                            legend.labs = comparisons,
                            legend.title = 'Mutation order',
                            legend = "none",
                            linetype = c("solid", "dashed", "solid"),
                            ggtheme = theme_cowplot())
    
    plot_list[[j]] = surv_plot
    
    print(surv_plot)
    png(filename = paste("~/Desktop/MetaAML_results/Figure_4/survival_by_muation_category_ordering/pngs/",g1,"_",g2, ".png", sep = ""), res = 300, width = 3.5, height = 3.5, units = "in")
    #
    surv_plot
    print(surv_plot)
    dev.off()
    
  }
}

temp_final_hr_categories_order = as.data.frame(do.call(rbind, results_list))
temp_final_hr_categories_order$fdr = p.adjust(temp_final_hr_categories_order$log_rank_p, method = "fdr")


# summary plot of all survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = FALSE,
                            ncol = 4, nrow = 3)
ggsave("~/Desktop/MetaAML_results/Figure_4/Supplimental/summary_survival_plot_grid.pdf", plots, width = 15, height = 10)



# forest plot ####
temp_final_hr_categories_order$categories = paste(temp_final_hr_categories_order$Category_1, temp_final_hr_categories_order$Category_2, sep = " -> ")
temp_final_hr_categories_order$sig_color = 0

for(i in 1:nrow(temp_final_hr_categories_order)){
  if(temp_final_hr_categories_order$log_rank_p[i] < 0.05){
    temp_final_hr_categories_order$sig_color[i] =1
  }
}

temp_final_hr_categories_order$sig_color = as.factor(temp_final_hr_categories_order$sig_color)

temp_final_hr_categories_order = subset(temp_final_hr_categories_order, temp_final_hr_categories_order$Upper_CI != "Inf")

temp_final_hr_categories_order$categories <- reorder(temp_final_hr_categories_order$categories, temp_final_hr_categories_order$HR)

temp_final_hr_categories_order$p_text = NA
for(i in 1:nrow(temp_final_hr_categories_order)){
  if(temp_final_hr_categories_order$log_rank_p[i] < 0.05){
    temp_final_hr_categories_order$p_text[i] = as.numeric(temp_final_hr_categories_order$log_rank_p[i])
  }
}
temp_final_hr_categories_order$p_text = round(temp_final_hr_categories_order$p_text, 3)
temp_final_hr_categories_order$fdr = round(temp_final_hr_categories_order$fdr, 1)

temp_final_hr_categories_order$p_text = paste(" p =", temp_final_hr_categories_order$p_text, "; ", "q =", temp_final_hr_categories_order$fdr)

for(i in 1:nrow(temp_final_hr_categories_order)){
  if(temp_final_hr_categories_order$log_rank_p[i] > 0.05){
    temp_final_hr_categories_order$p_text[i] = ""
  }
}

temp_final_hr_categories_order$HR = as.numeric(temp_final_hr_categories_order$HR)
temp_final_hr_categories_order$Lower_CI = as.numeric(temp_final_hr_categories_order$Lower_CI)
temp_final_hr_categories_order$Upper_CI = as.numeric(temp_final_hr_categories_order$Upper_CI)

temp_final_hr_categories_order$categories = reorder(temp_final_hr_categories_order$categories, temp_final_hr_categories_order$HR)


ggplot(temp_final_hr_categories_order, aes(x = reorder(categories, -HR), y = HR, label = p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=1.5, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=2, linetype="dashed", color = "lightgrey") +
  geom_text(aes(categories, Upper_CI), hjust = 0, nudge_y = 0.1) +
  geom_pointrange(size = .75, stat = "identity", shape = 19, 
                  aes(x = categories, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "#737373", "1" = "#762a83"))+
  theme_cowplot() +
  ylab("Hazard Ratio")+
  ylim(0,4.5) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/gene_category_ordering_hr_forest_plot_de_novo_5.pdf", dpi = 300, width = 6.5, height = 4.5, units = "in")




# scatterplot matrix for frequent pairwise genotypes ####

load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")

scatterplot_data = subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 150 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(scatterplot_data)){
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "ITD"){
    scatterplot_data$Gene[i] <- "FLT3-ITD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "SNV"){
    scatterplot_data$Gene[i] <- "FLT3-TKD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "Deletion"){
    scatterplot_data$Gene[i] <- "FLT3-TKD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "INDEL"){
    scatterplot_data$Gene[i] <- "FLT3-ITD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "other"){
    scatterplot_data$Gene[i] <- "FLT3-TKD"
  }
}


# select useful columns
scatterplot_data = select(scatterplot_data, Sample, Gene, VAF_CN_corrected)

# remove FLT3-ITD rows because they have no VAF data
# scatterplot_data = subset(scatterplot_data, Gene != "FLT3-ITD")

# function to make a scatterplot for each pairwise genes entered
vaf_scatterplot_function <- function(gene_1, gene_2){
  
  # select mutations of interest
  sub_gene_1 <- subset(scatterplot_data, scatterplot_data$Gene == gene_1)
  sub_gene_2 <- subset(scatterplot_data, scatterplot_data$Gene == gene_2)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  if(nrow(gene_1_and_2) >= 5){
    
    gene_1_and_2=as.data.frame(gene_1_and_2) %>% distinct(Sample, VAF_CN_corrected.x, VAF_CN_corrected.y, .keep_all = F)
    
    gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF_CN_corrected.x),]
    gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
    
    # add point color columns for visualizing clonal/subclonal trends
    gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF_CN_corrected.x - gene_1_and_2$VAF_CN_corrected.y)
    
    gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
    
    if(nrow(gene_1_and_2) >= 5){
      
      # define point color
      gene_1_and_2$Clonality <- NA
      gene_1_and_2$Clonality_color <- NA
      
      for(i in 1:nrow(gene_1_and_2)){
        if(gene_1_and_2$vaf_ratio[i] <= 5 & gene_1_and_2$vaf_ratio[i] >= -5){
          gene_1_and_2$Clonality[i] <- "Ambiguous"
          gene_1_and_2$Clonality_color[i] = "Ambiguous"
        }
        if(gene_1_and_2$vaf_ratio[i] > 5){
          gene_1_and_2$Clonality[i] <- paste(gene_1, "first", sep = " ")
          gene_1_and_2$Clonality_color[i] = "gene1_first"
        }
        if(gene_1_and_2$vaf_ratio[i] < -5){
          gene_1_and_2$Clonality[i] <- paste(gene_2, "first", sep = " ")
          gene_1_and_2$Clonality_color[i] = "gene2_first"
        }
      }
      
      gene_1_and_2$Clonality_color <- as.factor(gene_1_and_2$Clonality_color)
      
      
      # if(nrow(gene_1_and_2) > 4){
      # make the scatterplot
      scatter_plot = ggplot(gene_1_and_2, aes(x = gene_1_and_2$VAF_CN_corrected.y, y = gene_1_and_2$VAF_CN_corrected.x)) +
        geom_point(aes(color = Clonality_color), size =3, alpha = 0.75) +
        theme_cowplot() +
        scale_x_continuous(
          limits = c(0,100),
          labels = scales::number_format(accuracy = 10)) +
        scale_y_continuous(
          limits = c(0,100),
          labels = scales::number_format(accuracy = 10)) +
        xlab(paste(gene_2, " VAF", sep ="")) +
        ylab(paste(gene_1, " VAF", sep ="")) +  
        scale_color_manual(values = c("Ambiguous" = "#374E55FF", "gene2_first" = "#8c510a", "gene1_first" = "#01665e")) +
        geom_abline(intercept = 5, slope = (1), color="#969696",
                    linetype="dashed", size=.5)+
        geom_abline(intercept = -5, slope = (1), color="#969696",
                    linetype="dashed", size=.5)+
        geom_point(shape = 1, size =  3,colour = "black") +
        theme(legend.position = "none", 
              legend.title = element_blank()
        ) 
      # }    
      
      print(scatter_plot)
    }
  }
}

# find all pairwise gene comparisons
scatterplot_genes=as.data.frame(t(combn(as.vector(unique(scatterplot_data$Gene)),2)))

scatterplot_genes <- scatterplot_genes[order(scatterplot_genes$V1),]

# make list to populate plots into 
sc_plots = list()

for(j in 1:nrow(scatterplot_genes)){
  print(j)
  g1 = as.character(scatterplot_genes[j,1])
  g2 = as.character(scatterplot_genes[j,2])
  
  sc_plots[[j]] = vaf_scatterplot_function(gene_1 = g1, gene_2 = g2)
}

sc_plots = list.clean(sc_plots)

p = cowplot::plot_grid(plotlist = sc_plots, 
                       ncol = 8, align = "hv",   
                       rel_heights = c(1,1),
                       rel_widths = c(1,1))

ggsave(filename ="~/Desktop/MetaAML_results/Figure_4/Supplimental/VAF_scatterplot_matrix.png", dpi = 150, width = 17, height = 22, units = "in")

