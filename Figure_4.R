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

final_data_matrix = final_data_matrix_2

# pairwise ordering ####
# calculate the number of cases for co-occuring mutations
sub <- subset(final_data_matrix, final_data_matrix$mut_freq_gene >= 50 & final_data_matrix$Gene != "MLL" & final_data_matrix$Subset == "de_novo")

sub$Gene = as.character(sub$Gene)

# annotate FLT3 ITD and TKD respectively
for(i in 1:nrow(sub)){
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "SNV"){
    sub$Gene[i] <- "FLT3-TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] != "SNV" & sub$variant_type[i] != "Unknown"){
    sub$Gene[i] <- "FLT3-ITD"
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
  geom_point(aes(size = temp_dat_final_melted$number_1_and_2, color = fraction_1_then_2), shape = 15, stat = "identity") +
  geom_point(aes(size = temp_dat_final_melted$number_1_and_2), shape = 0, color = "#374E55FF") +
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
if (!require('BradleyTerryScalable')) install.packages('BradleyTerryScalable'); library('BradleyTerryScalable')
if (!require('Matrix.utils')) install.packages('Matrix.utils'); library('Matrix.utils')

# create the requred BT dataframe
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 75 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo" & final_data_matrix_2$mut_freq_pt > 1)

sub$Gene = as.character(sub$Gene)

# annotate FLT3 ITD and TKD respectively
for(i in 1:nrow(sub)){
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "SNV"){
    sub$Gene[i] <- "FLT3-TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] != "SNV" & sub$variant_type[i] != "Unknown"){
    sub$Gene[i] <- "FLT3-ITD"
  }
}
sub=subset(sub, sub$Gene != "FLT3")

genes <- as.data.frame(unique(sub$Gene))
results_list <- list()
n=1

for(i in 1:nrow(genes)){
  # print(i)
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
      
      gene_1_and_2$vaf_ratio <- as.numeric(gene_1_and_2$vaf_ratio)
      
      gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
      
      if(nrow(gene_1_and_2) > 0){
        # define point color
        gene_1_and_2$Clonality <- NA
        
        for(k in 1:nrow(gene_1_and_2)){
          if(gene_1_and_2$vaf_ratio[k] <= 0.05 & gene_1_and_2$vaf_ratio[k] >= -0.05){
            gene_1_and_2$Clonality[k] <- 0
          }
          if(gene_1_and_2$vaf_ratio[k] > 0.05){
            gene_1_and_2$Clonality[k] <- 1
          }
          if(gene_1_and_2$vaf_ratio[k] < -0.05){
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



TB_btdata <- btdata(temp_final, return_graph = T)
library(igraph)
par(mar = c(0, 0, 0, 0) + 0.1)  
plot.igraph(TB_btdata$graph, vertex.size = 28, edge.arrow.size = 0.5) 

summary(TB_btdata)

TB_btdata_fit <- btfit(TB_btdata, 1)
BT_MetaAML_mutation_ordering=as.data.frame(summary(TB_btdata_fit, SE = TRUE)$item_summary)

# add functional category to the mutations for visualization purposes
BT_MetaAML_mutation_ordering$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(BT_MetaAML_mutation_ordering)){
  if(BT_MetaAML_mutation_ordering$item[i] %in% DNA_methylation){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "DNA Methylation"
  }
  if(BT_MetaAML_mutation_ordering$item[i] %in% Chromatin_cohesin){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(BT_MetaAML_mutation_ordering$item[i] %in% RTK_RAS_Signaling){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(BT_MetaAML_mutation_ordering$item[i] %in% Splicing){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Splicing"
  }
  if(BT_MetaAML_mutation_ordering$item[i] %in% Transcription){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Transcription"
  }
  if(BT_MetaAML_mutation_ordering$item[i] == "NPM1"){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "NPM1"
  }
  if(BT_MetaAML_mutation_ordering$item[i] %in% Tumor_suppressors){
    BT_MetaAML_mutation_ordering$mutation_category[i] <- "Tumor suppressors"
  }
}


b=ggplot(BT_MetaAML_mutation_ordering,aes(reorder(factor(item), estimate),
                                          y=estimate,ymin=(estimate-SE),ymax=(estimate+SE))) +
  geom_pointrange(size = 0.75,
                  aes(x = reorder(factor(item), estimate), ymin = (estimate-SE), ymax = (estimate+SE), y = estimate, color = mutation_category)) +
  # geom_line(size = 1, linetype = 2,
  #           aes(x = reorder(factor(item), estimate), ymin = 1, ymax = estimate,
  #               color = mutation_category)) +
  scale_color_manual(values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF"), name = "Mutation Category") +
  scale_y_continuous(position = "right") +
  ylab("Point Estimate + 95% CI")+
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank()) +
  theme(legend.position = c(0.5, 0.85))   + coord_flip()+ scale_y_reverse()

print(b)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/bradley_terry_order_de_novo.pdf", width = 7.5, height = 6, units = "in")


#define the ordering based on the global pairwise analysis
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 75 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo" & final_data_matrix_2$mut_freq_pt > 1)

sub$Gene = as.character(sub$Gene)

# annotate FLT3 ITD and TKD respectively
for(i in 1:nrow(sub)){
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "SNV"){
    sub$Gene[i] <- "FLT3-TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] != "SNV" & sub$variant_type[i] != "Unknown"){
    sub$Gene[i] <- "FLT3-ITD"
  }
}
sub=subset(sub, sub$Gene != "FLT3")
gene_order <- as.list(BT_MetaAML_mutation_ordering$item)
gene_order <- rev(gene_order)
sub <- setDT(sub)[Gene %in% gene_order]
sub$Gene <- factor(sub$Gene, levels = gene_order)

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


c = ggplot(sub, aes(y = Gene,  x = VAF_male_x, fill = mutation_category, height = ..density..)) +
  geom_density_ridges(
    alpha = 1, stat = "density"
  ) +
  ylab(label = NULL) +
  xlab(label = "VAF") +
  xlim(0,1) +
  theme(legend.title = element_text()) +
  scale_fill_manual(name = "", values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF")) +
  theme(legend.position = "none")  

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/vaf_distribution_order_by_time.png", dpi = 300, width = 5, height = 5, units = "in")



ggarrange(a,c,b,
          ncol = 3, nrow = 1, widths = c(4,1.5, 2))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/pairwise_ordering_panels.png", dpi = 300, width = 17, height = 7.5, units = "in")


# pairwise scatterplot function ####
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
# final_data_matrix <- read.csv("~/Desktop/MetaAML_results/final_data_matrix.csv")


# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "ITD"){
    final_data_matrix$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "SNV"){
    final_data_matrix$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "Deletion"){
    final_data_matrix$Gene[i] <- "FLT3-TKD"
  }
  if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "INDEL"){
    final_data_matrix$Gene[i] <- "FLT3-ITD"
  }
}

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
    if(gene_1_and_2$vaf_ratio[i] <= .05 & gene_1_and_2$vaf_ratio[i] >= -.05){
      gene_1_and_2$Clonality[i] <- "Ambiguous"
    }
    if(gene_1_and_2$vaf_ratio[i] > .05){
      gene_1_and_2$Clonality[i] <- paste(gene_x, "first", sep = " ")
    }
    if(gene_1_and_2$vaf_ratio[i] < -.05){
      gene_1_and_2$Clonality[i] <- paste(gene_y, "first", sep = " ")
    }
  }
  
  if(nrow(gene_1_and_2) > 4){
    # make the scatterplot
    
    scatter_plot = ggplot(gene_1_and_2, aes(x = gene_1_and_2$VAF.y, y = gene_1_and_2$VAF.x)) +
      # geom_rect(mapping=aes(xmin=-Inf, xmax=.5, ymin=-Inf, ymax=.5), fill = "lightgrey", color=NA, alpha=.5) +
      geom_point(aes(color = Clonality), size = 3, alpha = 0.75) +
      xlim(0,(max(gene_1_and_2$VAF.y) + .05))+
      ylim(0,(max(gene_1_and_2$VAF.x) + .05))+
      xlab(paste(gene_y, " VAF", sep ="")) +
      ylab(paste(gene_x, " VAF", sep ="")) +  
      
      scale_color_manual(values = c("#374E55FF","#8c510a","#01665e")) +
      geom_abline(intercept = .05, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_abline(intercept = -.05, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_point(shape = 1, size =  3,colour = "black") +
      theme(plot.title = element_text(hjust = 0.5, paste(gene_y, "vs.", gene_x, sep = "")), 
            legend.position = "right", 
            legend.title = element_blank(), 
            legend.text = element_text(size=10), 
            axis.title.x = element_text(size=10), 
            axis.text.x = element_text(size=8),
            axis.title.y = element_text(size=10),
            axis.text.y = element_text(size=8)) 
    
    print(scatter_plot)
    
    if(save_plot == T){
      ggsave(filename =paste("~/Desktop/MetaAML_results/Figure_4/",gene_x, "_",gene_y,"_scatterplot.png", sep = ""), dpi = 300, width = 4.5, height = 3, units = "in")
    }
  }
}


vaf_scatterplot_function(pt_subset = "De novo", gene_1_2 = c("NRAS", "PTPN11"), save_plot = T)




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


load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_sub$Time_to_OS <- (final_data_matrix_sub$Time_to_OS/365)

# individual gene pairwise survival ####
# select recurrent mutations
genes = subset(final_data_matrix_sub, final_data_matrix_sub$mut_freq_gene >= 25)
genes=as.data.frame(t(combn(as.vector(unique(genes$Gene)),2)))

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
  gene_1_and_2 = unique(select(gene_1_and_2, Sample, Censor.x, Time_to_OS.x, VAF.x, VAF.y))
  
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF.x),]
  gene_1_and_2 = gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  if(nrow(gene_1_and_2) >= 10){
    # add point color columns for visualizing clonal/subclonal trends
    gene_1_and_2$vaf_difference <- (gene_1_and_2$VAF.x - gene_1_and_2$VAF.y)
    
    gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_difference), ] 
    
    # define point color
    gene_1_and_2$Clonality <- NA
    
    for(i in 1:nrow(gene_1_and_2)){
      if(gene_1_and_2$vaf_difference[i] <= 0.05 & gene_1_and_2$vaf_difference[i] >= -0.05){
        gene_1_and_2$Clonality[i] <- 0
      }
      if(gene_1_and_2$vaf_difference[i] > 0.05){
        gene_1_and_2$Clonality[i] <- 1
      }
      if(gene_1_and_2$vaf_difference[i] < -0.05){
        gene_1_and_2$Clonality[i] <- 2
      }
    }
    
    # minor clone for gene 1 vs clonal
    gene_1_and_clonal <- subset(gene_1_and_2, gene_1_and_2$Clonality == 2)
    n_dist1 <- n_distinct(gene_1_and_clonal$Clonality)
    
    gene_2_and_clonal <- subset(gene_1_and_2, gene_1_and_2$Clonality == 1)
    n_dist2 <- n_distinct(gene_2_and_clonal$Clonality)
    
    gene_1_and_2 <- subset(gene_1_and_2, gene_1_and_2$Clonality != 0)
    n_dist_1_2 <- n_distinct(gene_1_and_2$Clonality)
    
    # calculate the p value for surival between patients with subclonal mutations in gene 1 compared to those with subclonal mutations in gene 2
    if(nrow(gene_2_and_clonal) >= 5 & nrow(gene_1_and_clonal) > 5){
      print(i)
      gene_1_and_2$Censor.x = as.numeric(gene_1_and_2$Censor.x)
      
      # hazard ratio compiled table
      model <- coxph( Surv(Time_to_OS.x, Censor.x) ~ Clonality,
                      data = gene_1_and_2 )
      
      forest_plot_data=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                         factor_id_sep = ":", sort_by = NULL)
      # extract the log-rank p-value for the individual comparisons
      forest_plot_data$log_rank_p = summary(model)$sctest[3]
      
      forest_plot_data$gene_1 = gene_1
      forest_plot_data$gene_2 = gene_2
      
      results_list[[n]] <- forest_plot_data
      n=n+1   
    }
  }
}


temp_final_hr = as.data.frame(do.call(rbind, results_list))
temp_final_hr$q_value <- p.adjust(temp_final_hr$p, method = "fdr")


write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Figure_4/survival_by_vaf_ordering_co_occuring_mutations.csv", row.names = F)


temp_final_hr = subset(temp_final_hr, temp_final_hr$Upper_CI != "Inf")
temp_final_hr$HR = as.numeric(temp_final_hr$HR)
temp_final_hr$gene1_gene2 = as.character(paste(temp_final_hr$gene_2, temp_final_hr$gene_1, sep = " -> "))
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
temp_final_hr$p_text = paste("p =", temp_final_hr$p_text)

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] > 0.05){
    temp_final_hr$p_text[i] = ""
  }
}
temp_final_hr$sig_color = as.factor(temp_final_hr$sig_color)
temp_final_hr$gene1_gene2 = as.character(temp_final_hr$gene1_gene2)
temp_final_hr$gene1_gene2 = reorder(temp_final_hr$gene1_gene2, temp_final_hr$HR)

ggplot(temp_final_hr, aes(x = gene1_gene2, y = HR, label = temp_final_hr$p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(gene1_gene2, Upper_CI), hjust = 0, nudge_y = 1) +
  geom_pointrange(size = .75, stat = "identity", shape = 21, fill = "white",
                  aes(x = gene1_gene2, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  geom_errorbar(size = 1, aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color, width = .5)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "darkred"))+
  ylab("Hazard Ratio")+
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip() +
  theme_cowplot()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/gene_order_hr_forest_plot_de_novo_05.pdf", dpi = 300, width = 5, height = 8, units = "in")






# mutation categories ####
# because there are still too few cases on a pairwise basis to perform survival analysis, I will now look at pairwise occurence more broadly in terms of mutation categories
# sine it appears that epigenetic mutations occur first and proliferation hits last, I will find all patients who have mutations in these two categories, determine the order of aquizition, and performe survival analysis between patients with different order of aquizition
ddir("~/Desktop/MetaAML_results/Figure_4/survival_by_muation_category_ordering/pngs")

load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_sub$Time_to_OS <- (final_data_matrix_sub$Time_to_OS/365)


final_data_matrix_sub$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(final_data_matrix_sub)){
  if(final_data_matrix_sub$Gene[i] %in% DNA_methylation){
    final_data_matrix_sub$mutation_category[i] <- "DNA Methylation"
  }
  if(final_data_matrix_sub$Gene[i] %in% Chromatin_cohesin){
    final_data_matrix_sub$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(final_data_matrix_sub$Gene[i] %in% RTK_RAS_Signaling){
    final_data_matrix_sub$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(final_data_matrix_sub$Gene[i] %in% Splicing){
    final_data_matrix_sub$mutation_category[i] <- "Splicing"
  }
  if(final_data_matrix_sub$Gene[i] %in% Transcription){
    final_data_matrix_sub$mutation_category[i] <- "Transcription"
  }
  if(final_data_matrix_sub$Gene[i] == "NPM1"){
    final_data_matrix_sub$mutation_category[i] <- "NPM1"
  }
  if(final_data_matrix_sub$Gene[i] %in% Tumor_suppressors){
    final_data_matrix_sub$mutation_category[i] <- "Tumor suppressor"
  }
}


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
    
    sub = sub[order(sub$Gene, -sub$VAF_male_x),]
    sub = sub[!duplicated(sub$Gene),]
    sub = sub[!duplicated(sub$mutation_category),]
    
    vaf_cat_1=sub$VAF_male_x[1]
    vaf_cat_2=sub$VAF_male_x[2]
    
    if(!is.na(vaf_cat_1) & !is.na(vaf_cat_2)){
      
      diff=as.numeric(vaf_cat_1-vaf_cat_2)
      
      if(diff >= 0.05){
        o=1
      } 
      if(diff <= -0.05){
        o=3
      } 
      if(diff > -0.05 & diff < 0.05){
        o=2
      } 
      for(k in 1:nrow(cat_1_cat_2)){
        if(cat_1_cat_2$Sample[k] == sample){
          cat_1_cat_2$order[k]=o
        }
      }
    }
  }
  
  # create the survival data object
  final=as.data.frame(cat_1_cat_2) %>% distinct(Sample, Censor, Time_to_OS, order, .keep_all = F)
  
  final=na.omit(final)
  final$Censor=as.numeric(final$Censor)
  final$Time_to_OS=as.numeric(final$Time_to_OS)
  
  final$OS <- with(final, Surv(Time_to_OS, Censor == 1))
  
  OS <- survfit(OS ~ order, data = final, conf.type = "log-log")
  # hazard ratio compiled table
  model <- coxph( Surv(Time_to_OS, Censor) ~ order,
                  data = final )
  
  forest_plot_data=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                     factor_id_sep = ":", sort_by = NULL)
  forest_plot_data$category_1 = mut_cat[j,1]
  forest_plot_data$category_2 = mut_cat[j,2]
  
  # extract the log-rank p-value for the individual comparisons
  forest_plot_data$log_rank_p = summary(model)$sctest[3]
  
  results_list[[n]] <- forest_plot_data
  n=n+1 
  
  # replace space with underscore in the mutation categories in order to save the plot
  g1 <- gsub(" ", "_", g1)
  g1 <- gsub("/", "_", g1)
  g2 <- gsub(" ", "_", g2)
  g2 <- gsub("/", "_", g2)
  
  # specify colors to match those in other plots
  cat_colors = c("DNA Methylation first" = "#374E55FF", "Chromatin/Cohesin first" = "#DF8F44FF", "RTK/RAS Signaling first" = "#00A1D5FF", "Splicing first" = "#B24745FF", "Transcription first" = "#79AF97FF", "NPM1 first" = "#80796BFF", "Tumor suppressor first" = "#6A6599FF", "Ambiguous" = "lightgrey")
  
  # define the comparison labels for plotting
  comparisons = c(paste(mut_cat[j,2], "first"), "Ambiguous", paste(mut_cat[j,1], "first"))
  
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
                          pval.coord = c(0, 0),
                          conf.int = F,
                          censor = T,
                          surv.median.line = "none",
                          risk.table = F,
                          risk.table.title = "",
                          risk.table.fontsize = 4,
                          risk.table.height = .3,
                          risk.table.y.text = T,
                          break.time.by = 5,
                          risk.table.pos = c("out"),
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
                          legend = "right",
                          linetype = c("solid", "dashed", "solid"),
                          ggtheme = theme_cowplot())
  
  plot_list[[j]] = surv_plot
  
  print(surv_plot)
  png(filename = paste("~/Desktop/MetaAML_results/Figure_4/survival_by_muation_category_ordering/pngs/",g1,"_",g2, ".png", sep = ""), res = 300, width = 5, height = 3, units = "in")
  #
  surv_plot
  print(surv_plot)
  dev.off()
  
}

temp_final_hr_categories_order = as.data.frame(do.call(rbind, results_list))
temp_final_hr_categories_order$fdr = p.adjust(temp_final_hr_categories_order$p, method = "fdr")


# summary plot of all survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = TRUE,
                            ncol = 4, nrow = 6)
ggsave("~/Desktop/MetaAML_results/Data/Figure_4/Supplimental/summary_survival_plot_grid.pdf", plots, width = 10, height = 14)



# forest plot ####
temp_final_hr_categories_order$categories = paste(temp_final_hr_categories_order$category_1, temp_final_hr_categories_order$category_2, sep = " -> ")
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
    temp_final_hr_categories_order$p_text[i] = temp_final_hr_categories_order$log_rank_p[i]
  }
}
temp_final_hr_categories_order$p_text = round(temp_final_hr_categories_order$p_text, 3)
temp_final_hr_categories_order$p_text = paste("p =", temp_final_hr_categories_order$p_text)

for(i in 1:nrow(temp_final_hr_categories_order)){
  if(temp_final_hr_categories_order$log_rank_p[i] > 0.05){
    temp_final_hr_categories_order$p_text[i] = ""
  }
}

ggplot(temp_final_hr_categories_order, aes(x = reorder(categories, -HR), y = HR, label = temp_final_hr_categories_order$p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=1.5, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=2, linetype="dashed", color = "lightgrey") +
  geom_text(aes(categories, Upper_CI), hjust = 0, nudge_y = 0.1) +
  geom_pointrange(size = .75, stat = "identity", shape = 15, 
                  aes(x = categories, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "#737373", "1" = "#762a83"))+
  ylab("Hazard Ratio")+
  ylim(0,4) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip() +
  theme_cowplot()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_4/gene_category_ordering_hr_forest_plot_de_novo_5.pdf", dpi = 300, width = 9, height = 4.5, units = "in")

