# MetaAML_vaf_plots_and_by_mutation_order.R
#
# Brooks Benard
# bbenard@stanford.edu

if (!require('shiny')) install.packages('shiny'); library('shiny')
if (!require('scales')) install.packages('scales'); library('scales')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('UpSetR')) install.packages('UpSetR'); library('UpSetR')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('RCurl')) install.packages('RCurl'); library('RCurl')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
#

load("~/Desktop/MetaAML_results/final_data_matrix.RData")

co_occuring_clonality_function <- function(pt_subset, gene_1_2_3, order_gene, save_plots){
  
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
  gene_x <- gene_1_2_3[1]
  gene_y <- gene_1_2_3[2]
  gene_z <- gene_1_2_3[3]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$symbol == gene_x)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$symbol == gene_y)
  sub_gene_3 <- subset(final_data_matrix_sub, final_data_matrix_sub$symbol == gene_z)
  
  # select patients with both mutations
  gene_12 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  gene_123 <- inner_join(gene_12, sub_gene_3, by = "Sample")
  
  gene_123 <- na.omit(gene_123)
  gene_123 <- gene_123[!duplicated(gene_123[,c("Sample")]),]
  
  n_sample <- as.numeric(nrow(gene_123))
  
  gene_123_sub <- select(gene_123, Sample, VAF.x, VAF.y, VAF)
  colnames(gene_123_sub) <- c("Sample", paste(gene_x), paste(gene_y), paste(gene_z))
  
  
  gene_123_sub_m <- melt(gene_123_sub)
  
  gene <- order_gene
  
  gene_123_sub <- subset(gene_123_sub_m, gene_123_sub_m$variable == gene)
  gene_123_sub <- gene_123_sub[order(gene_123_sub$value, decreasing = T),]
  order <- as.data.frame(gene_123_sub$Sample)
  order$order_sample <- 1:nrow(order)
  colnames(order)[1] <- "Sample"
  
  gene_123_sub_m <- left_join(gene_123_sub_m, order, by = "Sample")
  
  three_gene_vaf_plot <-  ggplot(gene_123_sub_m, aes(x = reorder(Sample, order_sample), value)) +
    geom_point(aes(color = as.factor(variable)), size = 2.5, alpha = .75) +
    scale_color_manual(values = c("darkred", "darkblue", "darkgreen")) +
    ggtitle(paste("Clonal Patterns of Co-occuring\n",gene_x, ", ", gene_y, ", & ", gene_z, " Mutations", sep = "")) +
    xlab(paste("Samples (n = ", n_sample, ")", sep = "")) +
    ylab("VAF") +
    theme(axis.text.x=element_blank(),
                  legend.key = element_rect(fill='NA'),
                  legend.key.size = unit(0.4, 'cm'),
                  legend.title = element_blank(),
                  legend.position = "right",
                  legend.text = element_text(size=8, face="bold"))
  
  print(three_gene_vaf_plot)
    if(save_plots == T){
      
      ggsave(filename = paste("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/",gene_x, "_",gene_y,"_", gene_z,"mutation_order_by_vaf.png", sep = ""), dpi = 300, width = 7.5, height = 5, units = "in")
    }
}

co_occuring_clonality_function(pt_subset = "All", gene_1_2_3 = c("ASXL1", "RUNX1", "SRSF2"), order_gene = "ASXL1", save_plots = T)  


# if(n_samples > 15){
#   
# }
# 
#   # survival analysis
#   # define the mutation aquizition order
#   gene_123_sub$mut_order
#   for(i in 1:nrow(gene_123_sub)){
#     if(gene_123_sub$NPM1_VAF[i] > gene_123_sub$FLT3_VAF[i]){
#       gene_123_sub$mut_order[i] <- 1
#     }
#     if(gene_123_sub$FLT3_VAF[i] > gene_123_sub$NPM1_VAF[i]){
#       gene_123_sub$mut_order[i] <- 2
#     }
# }
# 
# 
#   for(i in 1:nrow(gene_123_sub)){
#     if(gene_123_sub$DNMT3A_VAF[i] > gene_123_sub$NPM1_VAF[i] & gene_123_sub$DNMT3A_VAF[i] > gene_123_sub$FLT3_VAF[i]){
#       gene_123_sub$mut_order[i] <- 1
#     }
#     if(gene_123_sub$NPM1_VAF[i] > gene_123_sub$FLT3_VAF[i]){
#       gene_123_sub$mut_order[i] <- 2
#     }
#     if(gene_123_sub$FLT3_VAF[i] > gene_123_sub$DNMT3A_VAF[i] & gene_123_sub$FLT3_VAF[i] > gene_123_sub$NPM1_VAF[i]){
#       gene_123_sub$mut_order[i] <- 3
#     }
#   }
# 
# 
# 
# 
# 
#   sub <- select(sequential_data, Sample, Censor,Time_to_OS)
#   gene_123_sub <- left_join(gene_123_sub, sub, by = "Sample")
# 
#   gene_123_sub <- unique(gene_123_sub)
# 
#   gene_123_sub$Time_to_OS <- (gene_123_sub$Time_to_OS/365)
# 
#   # create the survival data
#   gene_123_sub$OS <- with(gene_123_sub, Surv(Time_to_OS, Censor == 1))
# 
#   # create the survival objects used to plot kaplan-meyer curves
#   OS <- survfit(OS ~ mut_order, data = gene_123_sub, conf.type = "log-log")
# 
#   # plots the survival
#   surv_plot <- ggsurvplot(OS,
#                           data = gene_123_sub,
#                           log = (OS),
#                           log.rank.weights = c("survdiff"),
#                           pval = T,
#                           # test.for.trend = T,
#                           pval.method.size = 3,
#                           pval.coord = c(0, 0),
#                           conf.int = F,
#                           censor = F,
#                           surv.median.line = "none",
#                           risk.table = T,
#                           risk.table.title = "",
#                           risk.table.fontsize = 5,
#                           risk.table.height = .25,
#                           risk.table.y.text = T,
#                           break.time.by = 1,
#                           risk.table.pos = c("out"),
#                           palette = c("darkred", "darkblue", "darkgreen"),
#                           title = paste("Survival by ",gene_x, "_",gene_y,"_", gene_z,"Co-occurence", sep = ""),
#                           xlab = "Years",
#                           ylim = c(0, 1.0),
#                           ylab =  "Survival Probability",
#                           font.main = c(20, "plain", "black"),
#                           pval.size = 5,
#                           font.x = c(20, "plain", "black"),
#                           font.y =  c(20, "plain", "black"),
#                           font.legend = c(15, "plain"),
#                           font.tickslab = c(15, "plain", "black"),
#                           legend.labs = c("NPM1 first", "FLT3 first"),
#                           legend.title = "Mutation status",
#                           legend = "right",
#                           ggtheme = theme(plot.title = element_text(hjust = 0.5)))
# 
#   print(surv_plot)
# 
#   png(filename = paste("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/survival_by_mutation_order_",gene_x, "_",gene_y,"_", gene_z,".png", sep = ""), res = 300, width = 15, height = 10, units = "in")
# 
#   surv_plot
#   print(surv_plot)
#   dev.off()
#   }
#   }
# }

co_occuring_clonality_function(pt_subset = "All", gene_1_2_3 = c("ASXL1", "RUNX1", "SRSF2"), order_gene = "ASXL1", save_plots = T)
  
  
  
  
  
  
