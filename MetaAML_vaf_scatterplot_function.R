# MetaAML_vaf_scatterplot
#
# Brooks Benard
# bbenard@stanford.edu
# 02/03/2020
#

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')

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
              legend.position = "top", 
              legend.title = element_blank(), 
              legend.text = element_text(size=10), 
              axis.title.x = element_text(size=10), 
              axis.text.x = element_text(size=8),
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=8)) 

    print(scatter_plot)
    
    if(save_plot == T){
      ggsave(filename =paste("~/Desktop/MetaAML_results/Data/Figures/vaf_scatterplots/MetaAML_",gene_x, "_",gene_y,"_scatterplot.png", sep = ""), dpi = 300, width = 4, height = 3.5, units = "in")
    }
  }
}


vaf_scatterplot_function(pt_subset = "De novo", gene_1_2 = c("NRAS", "GATA2"), save_plot = T)
  