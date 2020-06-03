# MetaAML_vaf_scatterplot_survival
#
# Brooks Benard
# bbenard@stanford.edu
# 03/14/2020
#

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')

# load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")
final_data_matrix <- read.csv("~/Desktop/MetaAML_results/final_data_matrix.csv")

vaf_scatterplot_survival_function <- function(pt_subset, gene_1_2, ambiguous_curve){
  
  # subset to desired cohorts
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix
    label1 <- "All"
  }
  if(pt_subset == "De novo" | pt_subset == "De_novo" | pt_subset == "de_novo"){
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
  g1 <- gene_1_2[1]
  g2 <- gene_1_2[2]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == g1)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == g2)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  gene_1_and_2=as.data.frame(gene_1_and_2) %>% distinct(Sample, Censor.x, Time_to_OS.x, VAF.x, VAF.y, .keep_all = F)
  
  # make sure to only include one row per patient if there are two mutations in the same gene
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF.x),]
  gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  # add point color columns for visualizing clonal/subclonal trends
  gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF.x - gene_1_and_2$VAF.y)
  
  gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
  
  # define point color
  gene_1_and_2$Clonality <- NA
  
  # assign order
  for(i in 1:nrow(gene_1_and_2)){
    if(gene_1_and_2$vaf_ratio[i] <= .05 & gene_1_and_2$vaf_ratio[i] >= -.05){
      gene_1_and_2$Clonality[i] <- "Ambiguous"
    }
    if(gene_1_and_2$vaf_ratio[i] > .05){
      gene_1_and_2$Clonality[i] <- paste(g1, "first", sep = " ")
    }
    if(gene_1_and_2$vaf_ratio[i] < -.05){
      gene_1_and_2$Clonality[i] <- paste(g2, "first", sep = " ")
    }
  }
  
  if(ambiguous_curve == T){
    # p_val = F
    pal = c("#374E55FF","#8c510a","#01665e")
    lab = c("Ambiguous",paste(g2, "first",sep = " "), paste(g1, "first",sep = " "))
  }
  if(ambiguous_curve == F){
    gene_1_and_2 = subset(gene_1_and_2, gene_1_and_2$Clonality != "Ambiguous")
     # p_val = T
     pal = c("#8c510a","#01665e")
     lab = c(paste(g2, "first",sep = " "), paste(g1, "first",sep = " "))
  }

  # create the survival data object
  gene_1_and_2=na.omit(gene_1_and_2)
  gene_1_and_2$Time_to_OS.x <- (gene_1_and_2$Time_to_OS.x/365)
  gene_1_and_2$Time_to_OS.x=as.numeric(gene_1_and_2$Time_to_OS.x)
  
  # create survival object
  gene_1_and_2$OS <- with(gene_1_and_2, Surv(Time_to_OS.x, Censor.x == 1))
  OS <- survfit(OS ~ Clonality, data = gene_1_and_2, conf.type = "log-log")
  
  # find all the p-values for the different comparisons
  res <- pairwise_survdiff(Surv(Time_to_OS.x, Censor.x) ~ Clonality,
                           data = gene_1_and_2)
  p_val = round(res$p.value[2,2], 3)
  
  # plots the survival
  surv_plot <- ggsurvplot(OS,
                          data = gene_1_and_2,
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
                          palette = pal,
                          xlab = "Years",
                          ylim = c(0, 1.0),
                          ylab =  "Survival Probability",
                          font.main = c(15, "plain", "#252525"),
                          pval.size = 4,
                          font.x = c(12, "plain", "#252525"),
                          font.y =  c(12, "plain", "#252525"),
                          font.legend = c(12, "plain"),
                          font.tickslab = c(12, "plain", "#252525"),
                          legend.labs = lab,
                          legend.title = 'Mutation order',
                          legend = "right",
                          ggtheme = theme(plot.title = element_text(hjust = 0.5)))
  
    print(surv_plot)
    png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/vaf_scatterplots/survival_by_vaf_ordering/",g1,"_",g2,".png", sep = ""), res = 300, width = 4.5, height = 3, units = "in")
    surv_plot
    print(surv_plot)
    dev.off()
}
  
# run the function
vaf_scatterplot_survival_function(pt_subset = "De novo", gene_1_2 = c("NRAS", "GATA2"), ambiguous_curve = T)
  
  