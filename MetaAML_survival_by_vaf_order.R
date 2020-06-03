# Brooks Benard
# bbenard@stanford.edu
# 11.28.19
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


write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Data/survival_by_vaf_ordering/survival_by_vaf_ordering_co_occuring_mutations.csv", row.names = F)


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
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/gene_order_hr_forest_plot_de_novo_05.pdf", dpi = 300, width = 5, height = 8, units = "in")






# mutation categories ####
# because there are still too few cases on a pairwise basis to perform survival analysis, I will now look at pairwise occurence more broadly in terms of mutation categories
# sine it appears that epigenetic mutations occur first and proliferation hits last, I will find all patients who have mutations in these two categories, determine the order of aquizition, and performe survival analysis between patients with different order of aquizition
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
  # surv_plot <- ggsurvplot(OS,
  #                         data = final,
  #                         log = (OS),
  #                         log.rank.weights = c("survdiff"),
  #                         pval = p_val,
  #                         test.for.trend = F,
  #                         pval.method.size = 3,
  #                         pval.coord = c(0, 0.05),
  #                         conf.int = F,
  #                         censor = T,
  #                         surv.median.line = "none",
  #                         risk.table = F,
  #                         risk.table.title = "",
  #                         risk.table.fontsize = 4,
  #                         risk.table.height = .3,
  #                         risk.table.y.text = T,
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
  #                         ggtheme = theme(plot.title = element_text(hjust = 0.5)))
  # 
  # plot_list[[j]] = surv_plot
  # 
  # print(surv_plot)
  # png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_muation_category_ordering/pngs/",g1,"_",g2, ".png", sep = ""), res = 300, width = 5, height = 3, units = "in")
  # # 
  # surv_plot
  # print(surv_plot)
  # dev.off()

}

temp_final_hr_categories_order = as.data.frame(do.call(rbind, results_list))
temp_final_hr_categories_order$fdr = p.adjust(temp_final_hr_categories_order$p, method = "fdr")


# summary plot of all survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = TRUE,
                    ncol = 4, nrow = 6)
ggsave("~/Desktop/MetaAML_results/Data/Figures/survival_by_muation_category_ordering/summary_survival_plot_grid.pdf", plots, width = 10, height = 14)



# ### 
# final_data_matrix_sub
# 
# final_data_matrix_sub$Censor=as.numeric(final_data_matrix_sub$Censor)
# final_data_matrix_sub$Time_to_OS=as.numeric(final_data_matrix_sub$Time_to_OS)
# 
# final_data_matrix_sub$OS <- with(final_data_matrix_sub, Surv(Time_to_OS, Censor == 1))
# 
# OS <- survfit(OS ~ mutation_category, data = final_data_matrix_sub, conf.type = "log-log")
# 
# cat_colors = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressor" = "#6A6599FF")
# 
# 
# ggsurvplot(OS,
#            data = final_data_matrix_sub,
#            log = (OS),
#            log.rank.weights = c("survdiff"),
#            pval = F,
#            test.for.trend = F,
#            pval.method.size = 3,
#            pval.coord = c(0, 0.05),
#            conf.int = F,
#            censor = T,
#            surv.median.line = "none",
#            risk.table = F,
#            risk.table.title = "",
#            risk.table.fontsize = 4,
#            risk.table.height = .3,
#            risk.table.y.text = T,
#            break.time.by = 5,
#            risk.table.pos = c("out"),
#            palette = c("#DF8F44FF", "#374E55FF", "#80796BFF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF"),
#            xlab = "Years",
#            ylim = c(0, 1.0),
#            ylab =  "Survival Probability",
#            font.main = c(15, "plain", "#252525"),
#            pval.size = 4,
#            font.x = c(12, "plain", "#252525"),
#            font.y =  c(12, "plain", "#252525"),
#            font.legend = c(12, "plain"),
#            font.tickslab = c(12, "plain", "#252525"),
#            legend.title = 'Mutation category\noccuring first',
#            legend = "right",
#            ggtheme = theme(plot.title = element_text(hjust = 0.5)))













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
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/survival_by_muation_category_ordering/gene_category_ordering_hr_forest_plot_de_novo_5.pdf", dpi = 300, width = 9, height = 5, units = "in")
 

