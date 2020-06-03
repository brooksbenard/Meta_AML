## Brooks Benard
# bbenard@stanford.edu
# 11.29.19
# This script determines survival differences individual mutations in AML based on the patterns of clonality

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
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis'); library('survivalAnalysis')
library(maxstat)


load("~/Desktop/MetaAML_results/final_data_matrix.RData")
# save(final_data_matrix_2,  file = "~/Desktop/MetaAML_results/final_data_matrix_2.RData")
final_data_matrix_2_sub = subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 50 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")
final_data_matrix_2_sub$Time_to_OS <- (final_data_matrix_2_sub$Time_to_OS/365)

final_data_matrix_2_sub$Gene = as.character(final_data_matrix_2_sub$Gene)

for(i in 1:nrow(final_data_matrix_2_sub)){
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "SNV"){
    final_data_matrix_2_sub$Gene[i] = "FLT3-TKD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "ITD"){
    final_data_matrix_2_sub$Gene[i] = "FLT3-ITD"
  }
}

n=n_distinct(final_data_matrix_2_sub$Gene)
genes=data.frame(unique(final_data_matrix_2_sub$Gene))

# '%ni%' <- Negate('%in%')

results_list = list()
n=1
for(i in 1:nrow(genes)){
  gene=genes[i,1]
  mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
  mut_pts = distinct(mut_pts, Sample, Gene, VAF, VAF_male_x, Time_to_OS, Censor)
  mut_pts$hr_stratifier_vaf = 0
  mut_pts$hr_stratifier_vaf_text = "Late"
  
  mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
  mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
  
  for(j in 1:nrow(mut_pts)){
    if(!is.na(mut_pts$VAF_male_x[j])){
      if(mut_pts$VAF_male_x[j] >= 0.3){
        mut_pts$hr_stratifier_vaf[j] = 1
        mut_pts$hr_stratifier_vaf_text[j] = "Early"
        
      } 
    }
  }
  
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
                            risk.table = F,
                            risk.table.title = "",
                            risk.table.fontsize = 4,
                            risk.table.height = .3,
                            risk.table.y.text = T,
                            break.time.by = 5,
                            risk.table.pos = c("out"),
                            palette = c("Early" = "#6A6599FF", "Late" = "#80796BFF"),
                            xlab = "Years",
                            ylim = c(0, 1.0),
                            ylab =  "Survival Probability",
                            font.main = c(15, "plain", "#252525"),
                            pval.size = 4,
                            font.x = c(12, "plain", "#252525"),
                            font.y =  c(12, "plain", "#252525"),
                            font.legend = c(12, "plain"),
                            font.tickslab = c(12, "plain", "#252525"),
                            legend.labs = c("0" = "Early", "1" = "Late"),
                            legend.title = paste("Mutation\ntiming"),
                            legend = "right",
                            title = gene,
                            ggtheme = theme(plot.title = element_text(hjust = 0.5)))
    
    print(surv_plot)
    png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/",gene,"survival_by_VAF.png"), res = 300, width = 3.5, height = 3.5, units = "in")
    
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
    temp_final$sig_color[i] =1
  }
}

temp_final$sig_color = as.factor(temp_final$sig_color)

temp_final$fdr = p.adjust(temp_final$log_rank_p, method = "fdr")

temp_final$sig_color = as.factor(temp_final$sig_color)
# temp_final_hr_order_15$categories = as.factor(temp_final_hr_order_15$categories)
# temp_final_hr_order_15$HR = as.numeric(temp_final_hr_order_15$HR)

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
  geom_pointrange(size = 0.75, stat = "identity", shape = 21, fill = "white",
                  aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "darkred"))+
  ylab("Hazard Ratio")+
  ylim(0,7)+
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()


ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/gene_vaf_discrete_hr_forest_plot_de_novo_30.pdf", dpi = 300, width = 5, height = 7.5, units = "in")




  
  # optimal vaf cutoff ####
# now define the optimal VAF cutoff for each gene
# loop through all genes to calculate survival differences between major and minor VAFs for each gene

z <- 1
temp_list <- list()

results_list = list()
n=1

for(i in 1:nrow(genes)){
  print(i)
  gene <- genes[i,1]
  
  final_data_matrix_2_sub2 <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$Gene == gene)
  final_data_matrix_2_sub2 = na.omit(distinct(final_data_matrix_2_sub2, Sample, Gene, VAF, VAF_male_x, clonality, Time_to_OS, Censor))
  
  final_data_matrix_2_sub2 <- final_data_matrix_2_sub2[order(final_data_matrix_2_sub2$Sample, -final_data_matrix_2_sub2$VAF_male_x),]
  final_data_matrix_2_sub2= final_data_matrix_2_sub2[!duplicated(final_data_matrix_2_sub2$Sample),]
  
  if(n_distinct(final_data_matrix_2_sub2$clonality) > 1 & nrow(final_data_matrix_2_sub2) >= 15){
    
    # find the cutoff
    final_data_matrix_2_sub2$vaf_threshold <- NA
    
    mstat <- maxstat.test(Surv(Time_to_OS, Censor) ~ VAF_male_x, data=final_data_matrix_2_sub2, 
                          smethod="LogRank", pmethod="exactGauss", 
                          abseps=0.01)
    
    png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/", gene, "_vaf_logrank_plot.png"),res = 300, width = 5, height = 5, units = "in")
    plot(mstat)
    dev.off() 
    
    threshold = as.numeric(mstat$estimate)
    
    for(i in 1:nrow(final_data_matrix_2_sub2)){
      if(!is.na(final_data_matrix_2_sub2$VAF_male_x[i])){
        if(final_data_matrix_2_sub2$VAF_male_x[i] >= threshold){
          final_data_matrix_2_sub2$vaf_threshold[i] = "over"
        }
        if(final_data_matrix_2_sub2$VAF_male_x[i] < threshold){
          final_data_matrix_2_sub2$vaf_threshold[i] = "under"
        }
      }
    }
    
    if(n_distinct(final_data_matrix_2_sub2$vaf_threshold) > 1){
      # create the survival data object
      final_data_matrix_2_sub2$OS <- with(final_data_matrix_2_sub2, Surv(Time_to_OS, Censor ==1))
      
      OS <- survfit(OS ~ vaf_threshold, data = final_data_matrix_2_sub2, conf.type = "log-log")
      
      p <- surv_pvalue(OS)$pval
      
      # print(p)
      # store the results in a dataframe
      temp <- data.frame(matrix(NA, nrow = 1, ncol = 3))
      names(temp) <- c("Gene", "vaf_threshold", "p_value")
      
      temp[1,1] <- gene
      temp[1,2] <- threshold
      temp[1,3] <- p
      
      z <- z + 1
      temp_list[[z]] <- temp
      
      
      # find the different p-values for the different comparisons
      res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ vaf_threshold,
                               data = final_data_matrix_2_sub2)
      print(res)
      # cohorts <- c("Clonal", "Clonal LOH", "Subclonal", "Subclonal LOH")
      cohorts <- c("over", "under")
      
      title <- paste(threshold, "VAF\nthreshold",sep = " ")
      
      if(p <= 0.05){
        # plots the survival
        surv_plot <- ggsurvplot(OS,
                                data = final_data_matrix_2_sub2,
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
                                risk.table.y.text = T,
                                break.time.by = 5,
                                risk.table.pos = c("out"),
                                palette = c("over" = "#6A6599FF",  "under" = "#80796BFF"),
                                xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(15, "plain", "#252525"),
                                pval.size = 4,
                                font.x = c(12, "plain", "#252525"),
                                font.y =  c(12, "plain", "#252525"),
                                font.legend = c(12, "plain"),
                                font.tickslab = c(12, "plain", "#252525"),
                                legend.labs = cohorts,
                                legend.title = paste(title),
                                legend = "right",
                                title = gene,
                                ggtheme = theme(plot.title = element_text(hjust = 0.5)))
        
        print(surv_plot)
        png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/",gene,"survival_by_VAF.png"), res = 300, width = 5, height = 5, units = "in")
        
        surv_plot
        print(surv_plot)
        dev.off()
      }
      # summarize results from a Cox model
      final_data_matrix_2_sub2$vaf_threshold = ifelse(final_data_matrix_2_sub2$vaf_threshold == "over", 1,0)
      
      model <- coxph( Surv(Time_to_OS, Censor) ~ vaf_threshold,
                      data = final_data_matrix_2_sub2 )
      
      forest_plot_data=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                         factor_id_sep = ":", sort_by = NULL)
      forest_plot_data$gene = gene
      # extract the log-rank p-value for the individual comparisons
      forest_plot_data$log_rank_p = summary(model)$sctest[3]
      
      results_list[[n]] <- forest_plot_data
      n=n+1   
    }
  }
}

temp_final = as.data.frame(do.call(rbind, temp_list))
temp_final$q_value <- p.adjust(temp_final$q_value, method = "fdr")

write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/optimal_vaf_threshold_for_survival_prediction.csv", row.names = F)


# forrest plot of HRs ####
temp_final_hr = as.data.frame(do.call(rbind, results_list))

temp_final_hr$gene <- factor(temp_final_hr$gene, levels = temp_final_hr$gene[order(temp_final_hr$HR)])

temp_final_hr$sig_color = 0

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] < 0.05){
    temp_final_hr$sig_color[i] =1
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
temp_final_hr$p_text = paste("p =", temp_final_hr$p_text)

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$p[i] > 0.05){
    temp_final_hr$p_text[i] = ""
  }
}

ggplot(temp_final_hr, aes(x = reorder(gene, -HR), y = HR, label = temp_final_hr$p_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(gene, Upper_CI), hjust = 0, nudge_y = 0.5) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 21, fill = "white",
                  aes(x = gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "darkred"))+
  ylab("Hazard Ratio")+
  ylim(0,11.5)+
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()


ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/optimal_vaf_thresholds/gene_vaf_optimal_vaf_thresholds_hr_forest_plot_de_novo.pdf", dpi = 300, width = 5, height = 7.5, units = "in")

