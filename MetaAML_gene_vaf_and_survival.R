## MetaAML_gene_vaf_and_survival
## Brooks Benard
# bbenard@stanford.edu
# 04/08/2020
# This script looks at survival outcomes in the Meta-AML cohort based on 1) a static clonal/subclonal threshold, 2) a dynamic threshold optimized for each mutation, and 3) if there are survival trends based on VAF tertiles

# load required librarys
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
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis'); library('survivalAnalysis')
library(maxstat)
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')

# load the initial dataframe
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# subset to frequent mutations that are present in the de novo cohort
final_data_matrix_2_sub = subset(final_data_matrix, final_data_matrix$mut_freq_gene >= 50 & final_data_matrix$Gene != "MLL" & final_data_matrix$Subset == "de_novo")
final_data_matrix_2_sub$Time_to_OS <- (final_data_matrix_2_sub$Time_to_OS/365)
final_data_matrix_2_sub$Censor = as.numeric(final_data_matrix_2_sub$Censor)
final_data_matrix_2_sub$Gene = as.character(final_data_matrix_2_sub$Gene)

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

final_data_matrix_2_sub = distinct(final_data_matrix_2_sub, Sample, Gene, VAF, VAF_male_x, Time_to_OS, Censor)

n=n_distinct(final_data_matrix_2_sub$Gene)
genes=data.frame(unique(final_data_matrix_2_sub$Gene))


# static vaf threshold ####

results_list = list()
n=1
for(i in 1:nrow(genes)){
  gene=genes[i,1]
  mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
  
  mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
  mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
  
# set 30% as the VAF threshold for all mutations
  mut_pts$hr_stratifier_vaf = ifelse(mut_pts$VAF_male_x > 0.3, 1,0)
  mut_pts$hr_stratifier_vaf_text = ifelse(mut_pts$VAF_male_x > 0.3, "Clonal", "Subclonal")
  
  mut_pts = na.omit(mut_pts)
  
  # run the Cox model
  mut_pts$Time_to_OS=as.numeric(mut_pts$Time_to_OS)
  
  model <- coxph( Surv(Time_to_OS, Censor) ~ hr_stratifier_vaf,
                  data = mut_pts)
  
  forest_plot_data=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                     factor_id_sep = ":", sort_by = NULL)
  forest_plot_data$gene = genes[i,1]
  
  # extract the log-rank p-value for the individual comparisons
  forest_plot_data$log_rank_p = summary(model)$sctest[3]
  
  results_list[[n]] <- forest_plot_data
  n=n+1   
  
  p=forest_plot_data$log_rank_p
  
  # plot if significant
  if(p <= 0.05){
    # plots the survival
    mut_pts$OS <- with(mut_pts, Surv(Time_to_OS, Censor == 1))
    
    # mut_pts$hr_stratifier_vaf = as.numeric(mut_pts$hr_stratifier_vaf)
    
    OS <- survfit(OS ~ hr_stratifier_vaf_text, data = mut_pts, conf.type = "log-log")
    
    surv_plot <- ggsurvplot(OS,
                            data = mut_pts,
                            log = (OS),
                            log.rank.weights = c("survdiff"),
                            pval = T,
                            test.for.trend = F,
                            pval.method.size = 3,
                            pval.coord = c(0, 0),
                            conf.int = F,
                            censor = T,
                            surv.median.line = "none",
                            risk.table = T,
                            risk.table.title = "",
                            risk.table.fontsize = 4,
                            risk.table.height = .3,
                            risk.table.y.text = F,
                            break.time.by = 5,
                            risk.table.pos = c("out"),
                            palette = c("Clonal" = "#cb181d", "Subclonal" = "#3690c0"),
                            xlab = "Years",
                            ylim = c(0, 1.0),
                            ylab =  "Survival Probability",
                            font.main = c(15, "plain", "#252525"),
                            pval.size = 4,
                            font.x = c(12, "plain", "#252525"),
                            font.y =  c(12, "plain", "#252525"),
                            font.legend = c(12, "plain"),
                            font.tickslab = c(12, "plain", "#252525"),
                            legend.labs = c("Clonal" = "Clonal", "Subclonal" = "Subclonal"),
                            legend.title = paste("Mutation\nclonality"),
                            legend = "right",
                            title = gene,
                            ggtheme = theme(plot.title = element_text(hjust = 0.5)))
    
    print(surv_plot)
    png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/",gene,"survival_by_VAF.png"), res = 300, width = 5, height = 4.5, units = "in")
    
    surv_plot
    print(surv_plot)
    dev.off()
  }
}

temp_final = as.data.frame(do.call(rbind, results_list))

# temp_final = subset(temp_final, temp_final$gene != "U2AF1")

temp_final$gene <- factor(temp_final$gene, levels = temp_final$gene[order(temp_final$HR)])

temp_final$sig_color = 0

for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] < 0.05){
    temp_final$sig_color[i] = 1
  }
}

temp_final$sig_color = as.factor(temp_final$sig_color)

temp_final$fdr = p.adjust(temp_final$log_rank_p, method = "fdr")

temp_final$sig_color = as.factor(temp_final$sig_color)

temp_final = subset(temp_final, temp_final$Upper_CI != "Inf")

temp_final$categories <- reorder(temp_final$categories, temp_final$HR)

temp_final$p_text = NA
for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] < 0.05){
    temp_final$p_text[i] = temp_final$log_rank_p[i]
  }
}
temp_final$q_text = NA
for(i in 1:nrow(temp_final)){
  if(temp_final$fdr[i] < 0.2){
    temp_final$q_text[i] = temp_final$fdr[i]
  }
}
temp_final$p_text = round(temp_final$p_text, 3)
temp_final$q_text = round(temp_final$q_text, 2)

temp_final$p_q_text = paste("p =", temp_final$p_text, "; q =", temp_final$q_text)
temp_final$p_text = paste("p =", temp_final$p_text)

for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] > 0.05){
    temp_final$p_q_text[i] = ""
  }
  if(temp_final$log_rank_p[i] > 0.05){
    temp_final$p_text[i] = ""
  }
}

ggplot(temp_final, aes(x = reorder(gene, -HR), y = HR, label = temp_final$p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=3, linetype="dashed", color = "lightgrey") +
  geom_text(aes(gene, Upper_CI), hjust = 0, nudge_y = 0.1) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 15, 
                  aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  geom_errorbar(size = 1, aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color, width = .5)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "#762a83"))+
  ylab("Hazard Ratio")+
  # ylim(0,7)+
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()


ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/gene_vaf_discrete_hr_forest_plot_de_novo_30.pdf", dpi = 300, width = 5, height = 7, units = "in")




  
  # optimal vaf cutoff ####
# now define the optimal VAF cutoff for each gene
# loop through all genes to calculate survival differences between major and minor VAFs for each gene

z <- 1
temp_list <- list()

n=1
results_list = list()

# plot the optimal cutpoint for all genes
for(i in 1:nrow(genes)){
  print(i)
  gene <- as.character(genes[i,1])
  
  mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
  
  # make sure to remove rows with duplicate 
  mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
  mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
  
  if(nrow(mut_pts) >= 15){
  
    cut.point <- surv_cutpoint(
      mut_pts,
      time = "Time_to_OS",
      event = "Censor",
      variables = c("VAF_male_x")
    )

    p = plot(cut.point, palette = "aaas")
    
    print(p)
    png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/thresholds/",gene,"_optimal_VAF_threshold.png", sep = ""), res = 300, width = 5, height = 5, units = "in")
    
    p
    print(p)
    dev.off()
  }
}
    

# use the optimal cutpoint to plot the survival curves
for(i in 1:nrow(genes)){
      print(i)
      gene <- as.character(genes[i,1])
      
      mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
      
      # make sure to remove rows with duplicate 
      mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
      mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
      
      if(nrow(mut_pts) >= 15){
        # find the optimal cutpoint
        mut_pts$vaf_threshold <- NA
        
        mstat <- maxstat.test(Surv(Time_to_OS, Censor) ~ VAF_male_x, data=mut_pts,
                              smethod="LogRank", pmethod="exactGauss",
                              abseps=0.05)
        
        # png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/thresholds/", gene, "_vaf_logrank_plot.png",sep = ""),res = 300, width = 5, height = 5, units = "in")
        # plot(mstat)
        # dev.off()
        
        # extract the threshold value for future use
        threshold = as.numeric(mstat$estimate)
        
        mut_pts$vaf_threshold = ifelse(mut_pts$VAF_male_x > threshold, 1,0)
        mut_pts$vaf_threshold_text = ifelse(mut_pts$VAF_male_x > threshold, "Above", "Below")
        
        # remove incomplete data
        mut_pts = na.omit(mut_pts)
        
    # summarize results from a Cox model
    model <- coxph( Surv(Time_to_OS, Censor) ~ vaf_threshold,
                    data = mut_pts )
    
    forest_plot_data = cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                       factor_id_sep = ":", sort_by = NULL)
    forest_plot_data$gene = gene
    # extract the log-rank p-value for the individual comparisons
    forest_plot_data$log_rank_p = summary(model)$sctest[3]
    
    # append results into a list
    results_list[[n]] <- forest_plot_data
    n=n+1   
    
      # make the survival plot for significant cases
        if(n_distinct(mut_pts$vaf_threshold) > 1){
          # create the survival data object
          mut_pts$OS <- with(mut_pts, Surv(Time_to_OS, Censor ==1))
          
          mut_pts$vaf_threshold = as.numeric(mut_pts$vaf_threshold)
          
          OS <- survfit(OS ~ vaf_threshold_text, data = mut_pts, conf.type = "log-log")
          
          p <- surv_pvalue(OS)$pval

          # store the results in a dataframe
          temp <- data.frame(matrix(NA, nrow = 1, ncol = 3))
          names(temp) <- c("Gene", "vaf_threshold", "p_value")
          
          temp[1,1] <- gene
          temp[1,2] <- threshold
          temp[1,3] <- p
          
          z <- z + 1
          temp_list[[z]] <- temp
          
          thresh <- paste(threshold, "VAF\nthreshold",sep = " ")
      
              if(p <= 0.05){
                # plots the survival
                surv_plot <- ggsurvplot(OS,
                                        data = mut_pts,
                                        log = (OS),
                                        log.rank.weights = c("survdiff"),
                                        pval = T,
                                        test.for.trend = F,
                                        pval.method.size = 3,
                                        pval.coord = c(0, 0),
                                        conf.int = F,
                                        censor = T,
                                        surv.median.line = "none",
                                        risk.table = T,
                                        risk.table.title = "",
                                        risk.table.fontsize = 4,
                                        risk.table.height = .3,
                                        risk.table.y.text = F,
                                        break.time.by = 5,
                                        risk.table.pos = c("out"),
                                        palette = c("Above" = "#cb181d", "Below" = "#3690c0"),
                                        xlab = "Years",
                                        ylim = c(0, 1.0),
                                        ylab =  "Survival Probability",
                                        font.main = c(15, "plain", "#252525"),
                                        pval.size = 4,
                                        font.x = c(12, "plain", "#252525"),
                                        font.y =  c(12, "plain", "#252525"),
                                        font.legend = c(12, "plain"),
                                        font.tickslab = c(12, "plain", "#252525"),
                                        legend.labs = c("Above" = "Above", "Below" = "Below"),
                                        legend.title = paste(thresh),
                                        legend = "right",
                                        title = gene,
                                        ggtheme = theme(plot.title = element_text(hjust = 0.5)))
                
                print(surv_plot)
                png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/survival/",gene,"_survival_by_VAF.png", sep = ""), res = 300, width = 5, height = 5, units = "in")
                
                surv_plot
                print(surv_plot)
                dev.off()
      }
    }
  }
}

temp_final = as.data.frame(do.call(rbind, temp_list))
temp_final$q_value <- p.adjust(temp_final$p_value, method = "fdr")

write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/optimal_vaf_threshold_for_survival_prediction.csv", row.names = F)


# forrest plot of HRs ####
temp_final_hr = as.data.frame(do.call(rbind, results_list))

temp_final_hr$gene <- factor(temp_final_hr$gene, levels = temp_final_hr$gene[order(temp_final_hr$HR)])

temp_final_hr$sig_color = 0

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] < 0.05 & temp_final_hr$HR[i] < 1){
    temp_final_hr$sig_color[i] = 1
  }
  if(temp_final_hr$log_rank_p[i] < 0.05 & temp_final_hr$HR[i] > 1){
    temp_final_hr$sig_color[i] = 2
  }
}

temp_final_hr$sig_color = as.factor(temp_final_hr$sig_color)

temp_final_hr$fdr = p.adjust(temp_final_hr$log_rank_p, method = "fdr")

temp_final_hr$sig_color = as.factor(temp_final_hr$sig_color)

temp_final_hr = subset(temp_final_hr, temp_final_hr$Upper_CI != "Inf")

temp_final_hr$categories <- reorder(temp_final_hr$categories, temp_final_hr$HR)

temp_final_hr$p_text = NA
for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] < 0.05){
    temp_final_hr$p_text[i] = temp_final_hr$log_rank_p[i]
  }
}

temp_final_hr$p_text = round(temp_final_hr$p_text, 3)
temp_final_hr$p_text = ifelse(temp_final_hr$log_rank_p < 0.001, "p < 0.001",paste("p =", temp_final_hr$p_text))

# temp_final_hr$p_text = paste("p =", temp_final_hr$p_text)

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$p[i] > 0.05){
    temp_final_hr$p_text[i] = ""
  }
}

ggplot(temp_final_hr, aes(x = reorder(gene, -HR), y = HR, label = temp_final_hr$p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2.5, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=5, linetype="dashed", color = "lightgrey") +
  geom_text(aes(gene, Upper_CI), hjust = 0, nudge_y = 0.5, size = 3) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 15, 
                  aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  geom_errorbar(size = 1, aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color, width = .5)) +
  scale_color_manual(values = c("0" = "#737373", "1" = "#1b7837", "2" = "#762a83"))+
  ylab("Hazard Ratio")+
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 


ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/gene_vaf_optimal_vaf_thresholds_hr_forest_plot_de_novo.pdf", dpi = 300, width = 5, height = 6, units = "in")





# VAF tertiles ####
n=1
results_list = list()

for(i in 1:nrow(genes)){
  print(i)
  gene = as.character(genes[i,1])
  
  if(gene != "FLT3-ITD"){
    mut_pts=subset(final_data_matrix_2_sub, Gene == gene)
    
    mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
    mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
    
    # split the patients into vaf tertiles 
    mut_pts$tertile = ntile(mut_pts$VAF_male_x, 3)
    
    mut_pts = na.omit(mut_pts)
    
    
    # summarize results from a Cox model
    model <- coxph( Surv(Time_to_OS, Censor) ~ VAF_male_x,
                    data = mut_pts )
    
    forest_plot_data = cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                         factor_id_sep = ":", sort_by = NULL)
    forest_plot_data$gene = gene
    # extract the log-rank p-value for the individual comparisons
    forest_plot_data$log_rank_p = summary(model)$sctest[3]
    
    # append results into a list
    results_list[[n]] <- forest_plot_data
    n=n+1   
    
    
    
    mut_pts$OS <- with(mut_pts, Surv(Time_to_OS, Censor == 1))
    
    OS <- survfit(OS ~ tertile, data = mut_pts, conf.type = "log-log")
    
    surv_plot <- ggsurvplot(OS,
                            data = mut_pts,
                            log = (OS),
                            log.rank.weights = c("survdiff"),
                            pval = T,
                            test.for.trend = T,
                            pval.method.size = 3,
                            pval.coord = c(0, 0),
                            conf.int = F,
                            censor = T,
                            surv.median.line = "none",
                            risk.table = T,
                            risk.table.title = "",
                            risk.table.fontsize = 4,
                            risk.table.height = .3,
                            risk.table.y.text = F,
                            break.time.by = 5,
                            risk.table.pos = c("out"),
                            palette = c("#b35806", "lightgrey", "#542788"),
                            xlab = "Years",
                            ylim = c(0, 1.0),
                            ylab =  "Survival Probability",
                            font.main = c(15, "plain", "#252525"),
                            pval.size = 4,
                            font.x = c(12, "plain", "#252525"),
                            font.y =  c(12, "plain", "#252525"),
                            font.legend = c(12, "plain"),
                            font.tickslab = c(12, "plain", "#252525"),
                            legend.labs = c("1" = "Upper", "2" = "Middle", "3" = "Lower"),
                            legend.title = paste("VAF tertile"),
                            legend = "right",
                            title = gene,
                            ggtheme = theme(plot.title = element_text(hjust = 0.5)))
    
    print(surv_plot)
    
    png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/tertile_threshold/",gene,"survival_by_VAF_tertile.png"), res = 300, width = 5, height = 5.5, units = "in")
    
    surv_plot
    print(surv_plot)
    dev.off() 
  }
}


temp_final = as.data.frame(do.call(rbind, results_list))
temp_final$q_value <- p.adjust(temp_final$q_value, method = "fdr")

write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/optimal_vaf_threshold_for_survival_prediction.csv", row.names = F)


