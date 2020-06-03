
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('gganimate')) install.packages('gganimate'); library('gganimate')

load("~/Desktop/MetaAML_results/final_data_matrix.RData")
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



results_list = list()
n=1



for(i in 1:nrow(genes)){
  
  saveGIF({
  i = 1
  gene=genes[i,1]
  mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
  
  mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
  mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
  
  # vaf_range = range(mut_pts$VAF_male_x, na.rm = T)
  
  for(j in 2:19){
    cutpoint = j*.05
  
      mut_pts$threshold = cutpoint
      mut_pts$hr_stratifier_vaf = ifelse(mut_pts$VAF_male_x > cutpoint, 1,0)
      mut_pts$hr_stratifier_vaf_text = ifelse(mut_pts$VAF_male_x > cutpoint, "Early", "Late")
      
      mut_pts_final = na.omit(mut_pts)
      
      results_list[[n]] <- mut_pts_final
      n=n+1 
      
      }
  
  temp_final = as.data.frame(do.call(rbind, results_list))
  
  
    # plots the survival
  temp_final$OS <- with(temp_final, Surv(Time_to_OS, Censor == 1))
    
    for(k in 2:19){
      
      cutpoint = k*.05
      
      temp_final_sub = subset(temp_final, temp_final$threshold == cutpoint)
      
      temp_final_sub = distinct(temp_final_sub)
      
      OS <- survfit(OS ~ hr_stratifier_vaf_text, data = temp_final_sub, conf.type = "log-log")
      
      surv_plot <- ggsurvplot(OS,
                              data = temp_final_sub,
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
                              risk.table.y.text = F,
                              break.time.by = 5,
                              risk.table.pos = c("out"),
                              palette = c("Early" = "#b35806", "Late" = "#542788"),
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
                              legend.title = paste("VAF threshold\n", cutpoint),
                              legend = "right",
                              title = gene,
                              ggtheme = theme(plot.title = element_text(hjust = 0.5)))
      

      print(surv_plot)
      
      gene = gsub(" ", "", gene, fixed = TRUE)
      cutpoint = gsub(" ", "", cutpoint, fixed = TRUE)
      
      png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/gifs/",gene,"_survival_by_",cutpoint,"_VAF_threshold.png", sep = ""), res = 300, width = 5, height = 4.5, units = "in")
      
      surv_plot
      print(surv_plot)
      dev.off()
    }
}, interval = .2, movie.name="~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/gifs/test.gif")

if (!require('magick')) install.packages('magick'); library('magick')


png.files <- list.files("~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/gifs/", all.files = T) #Mention the number of files to read
GIF.convert <- function(x, output = "~/Desktop/MetaAML_results/Data/Figures/survival_by_vaf_per_gene/gifs/NRAS_animation.gif")#Create a function to read, animate and convert the files to gif
{
  image_read(x) %>%
    image_animate(fps = 1) %>%
    image_write(output)
}

GIF.convert(png.files)
