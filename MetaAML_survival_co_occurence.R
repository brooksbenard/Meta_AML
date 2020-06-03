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
    # dat_1 = subset(temp_sub_final, STRATA == 0 | STRATA == 1 | STRATA == 6)
    # dat_1$STRATA[dat_1$STRATA==6] = 1
    # 
    # dat_2 = subset(temp_sub_final, STRATA == 0 | STRATA == 2 | STRATA == 6)
    # dat_2$STRATA[dat_2$STRATA==6] = 2
    
    dat_3 = subset(temp_sub_final, STRATA == 1 | STRATA == 6)
    dat_4 = subset(temp_sub_final, STRATA == 2 | STRATA == 6)
    dat_5 = subset(temp_sub_final, STRATA == 0 | STRATA == 6)
    
    # run the Cox models
    # gene 1 vs wt
    # model_1 <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
    #                   data = dat_1)
    # forest_1=cox_as_data_frame(coxphsummary = model_1, unmangle_dict = NULL,
    #                            factor_id_sep = ":", sort_by = NULL)
    # forest_1$genes = paste(gene1)
    # forest_1$gene_tested = paste(gene1, "vs_wt", sep = "_")
    # 
    # # extract the log-rank p-value for the individual comparisons
    # forest_1$log_rank_p = summary(model_1)$sctest[3]
    # 
    # # gene 2 vs wt
    #       model_2 <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
    #                         data = dat_2)
    #       forest_2=cox_as_data_frame(coxphsummary = model_2, unmangle_dict = NULL,
    #                                  factor_id_sep = ":", sort_by = NULL)
    #       forest_2$genes = paste(gene2)
    #       forest_2$gene_tested = paste(gene2, "vs_wt", sep = "_")
          
          # extract the log-rank p-value for the individual comparisons
          # forest_2$log_rank_p = summary(model_2)$sctest[3]
          
    # co-mutated vs gene 1
     #              model_3 <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
     #                                data = dat_3)
     #              forest_3=cox_as_data_frame(coxphsummary = model_3, unmangle_dict = NULL,
     #                                         factor_id_sep = ":", sort_by = NULL)
     #              forest_3$gene_1 = gene1
     #              forest_3$gene_2 = gene2
     #              forest_3$gene_tested = paste(gene1, "vs_both", sep = "_")
     #              
     #              # extract the log-rank p-value for the individual comparisons
     #              forest_3$log_rank_p = summary(model_3)$sctest[3]
     #              
     # # co-mutated vs gene 2
     #                      model_4 <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
     #                                        data = dat_4)
     #                      forest_4=cox_as_data_frame(coxphsummary = model_4, unmangle_dict = NULL,
     #                                                 factor_id_sep = ":", sort_by = NULL)
     #                      forest_4$gene_1 = gene1
     #                      forest_4$gene_2 = gene2
     #                      forest_4$gene_tested = paste(gene2, "vs_both", sep = "_")
     #                      
     #                      # extract the log-rank p-value for the individual comparisons
     #                      forest_4$log_rank_p = summary(model_4)$sctest[3]
  
    # co-mutated vs wt
                          model_5 <- coxph( Surv(Time_to_OS, Censor) ~ STRATA,
                                            data = dat_5)
                          forest_5=cox_as_data_frame(coxphsummary = model_5, unmangle_dict = NULL,
                                                     factor_id_sep = ":", sort_by = NULL)
                          forest_5$gene_1 = gene1
                          forest_5$gene_2 = gene2
                          forest_5$gene_tested = paste(gene1, gene2, "vs_wt", sep = "_")
                          
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
                              ggtheme = theme(plot.title = element_text(hjust = 0.5)))

      plot_list[[i]] = surv_plot
      
      print(surv_plot)
    # 
    #   png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_co_occurence/",gene1, "_", gene2, "_co_occurence_survival.png", sep = ""), res = 300, width = 6, height = 4.5, units = "in")
    # 
    #   surv_plot
    #   print(surv_plot)
    #   dev.off()
    # }
    }
  }
}

# bind all data
temp_final_hr = as.data.frame(do.call(rbind, results_list))

# correct for mulitple hypothesis testing
temp_final_hr$q_value <- p.adjust(temp_final_hr$log_rank_p, method = "fdr")

# summary plot of all pairwise survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = TRUE,
                            ncol = 5, nrow = 8)
ggsave("~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/pairwise_survival_grid.pdf", plots, width = 17, height = 20)



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
ggsave("~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/indiviaual_mutations/survival_grid.pdf", plots, width = 16, height = 17)

  
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

pdf(file = "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/pairwise_survival_correlation_de_novo.pdf", width = 7.5, height = 7.5)

corrplot(freq_genes_matrix_hr,
         is.corr = F, 
         type="upper", 
         order="hclust",
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

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/indiviaual_mutations/forrest_plot_de_novo.pdf", dpi = 300, width = 7.5, height = 7.5, units = "in")



# scattterplot of HR and odds ratio ####
# read in the results from MetaAML_co_occurence_mutual_exlclusive_analysis.R
odds_ratio = read.csv("~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/oods_ratio_and_fishers_results.csv")
odds_ratio = select(odds_ratio, gene1, gene2, odds_ratio, fishers_q)
colnames(odds_ratio) = c("gene_1", "gene_2", "odds_ratio", "q_value_odds")

hazard_ratio = pairwise_data[,c(1:3, 5)]
colnames(hazard_ratio)[4] = "q_value_HR"

hr_odds = inner_join(odds_ratio, hazard_ratio, by = c("gene_1", "gene_2"))

# hr_odds$color = as.factor(ifelse(hr_odds$q_value_HR < 0.1, "HR q < 0.1","HR q > 0.1"))

hr_odds$color = "NS"

for(i in 1:nrow(hr_odds)){
  if(hr_odds$q_value_HR[i] < 0.1 & hr_odds$HR[i] > 1){
    hr_odds$color[i] = "HR > 1 & q < 0.1"
  }
  if(hr_odds$q_value_HR[i] < 0.1 & hr_odds$HR[i] < 1){
    hr_odds$color[i] = "HR < 1 & q < 0.1"
  }
}


ggplot(hr_odds, aes(x = log(odds_ratio), y = log(HR))) +
  geom_point(aes(color = color, shape = ), size = 3, alpha = 0.75) +
  xlab("log(Odds Ratio)") +
  ylab("log(Hazard Ratio)") +
  scale_color_manual(values = c("#1b7837",  "#762a83", "grey")) +
  geom_point(shape = 1, size =  3,colour = "black") +
  geom_hline(yintercept=log(1), linetype="dashed", color = "#d9d9d9") +
  geom_vline(xintercept=log(1), linetype="dashed", color = "#d9d9d9") +
  theme(
        legend.title = element_blank(), 
        legend.text = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(size=8),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=8)) 


ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/mutation_co_occurence/correlation_of_odds_ratio_and_HR_de_novo.pdf", dpi = 300, width = 5, height = 4, units = "in")
