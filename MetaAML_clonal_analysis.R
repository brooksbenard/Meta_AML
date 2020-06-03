# MetaAML_clonality_analysis
#
# Brooks Benard
# bbenard@stanford.edu
# 10/09/2019
#

### This script uses variant allele frequencies from MetaAML to understand clonal patterns at both genetic and survival levels

#### load required packages and data ####
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


# load the MetaAML results file 
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_2 <- read.csv("~/Desktop/MetaAML_results/final_data_matrix.csv")

# Assign clonality  by VAF ####
# assign simple clonal or subclonal calls to mutations based on the vaf and correct for x-linked genes in males
# filter vafs
final_data_matrix_2 = subset(final_data_matrix_2, final_data_matrix_2$VAF_male_x >= 0.01)

final_data_matrix_2$clonality <- NA

for(i in 1:nrow(final_data_matrix_2)){
  if(!is.na(final_data_matrix_2$VAF_male_x[i])){
    if(final_data_matrix_2$VAF_male_x[i] >=0.3){
      final_data_matrix_2$clonality[i] <- "Major"
    }
    if(final_data_matrix_2$VAF_male_x[i] < 0.3 & final_data_matrix_2$VAF_male_x[i] > 0){
      final_data_matrix_2$clonality[i] <- "Minor"
    }  
  }
}

# try to assign gradient vaf bins
# this is very arbitrary and I should use more rigorous methods (e.g. PyClones, SciClone, etc.) with the BeatAML and TCGA data in a future analysis
final_data_matrix_2$clonality_bin <- NA

for(i in 1:nrow(final_data_matrix_2)){
  if(!is.na(final_data_matrix_2$VAF_male_x[i])){
    # if(final_data_matrix_2$VAF_male_x[i] < 0.1 ){
    #   final_data_matrix_2$clonality_bin[i] <- "< 0.1"
    # }
    if(final_data_matrix_2$VAF_male_x[i] < 0.2){
      final_data_matrix_2$clonality_bin[i] <- "< 0.2"
    }
    if(final_data_matrix_2$VAF_male_x[i] >= 0.2 & final_data_matrix_2$VAF_male_x[i] < 0.4){
      final_data_matrix_2$clonality_bin[i] <- "0.2-0.4"
    }
    if(final_data_matrix_2$VAF_male_x[i] > 0.4 & final_data_matrix_2$VAF_male_x[i] <= 0.6){
      final_data_matrix_2$clonality_bin[i] <- "0.4-0.6"
    }
    # if(final_data_matrix_2$VAF_male_x[i] >= 0.4 & final_data_matrix_2$VAF_male_x[i] <= 0.5){
    #   final_data_matrix_2$clonality_bin[i] <- "0.4-0.5"
    # }
    # if(final_data_matrix_2$VAF_male_x[i] >= 0.5 & final_data_matrix_2$VAF_male_x[i] <= 0.6){
    #   final_data_matrix_2$clonality_bin[i] <- "0.5-0.6"
    # }
    if(final_data_matrix_2$VAF_male_x[i] > 0.6){
      final_data_matrix_2$clonality_bin[i] <- "> 0.6"
    }
  }
}


# count the number of "clones" per patient by finding the number of unique vaf bins
pts <- as.data.frame(unique(final_data_matrix_2$Sample))

z <- 1
temp_list <- list()

for(i in 1:nrow(pts)){
  pt <- as.character(pts[i,1])
  sub_dat <- subset(final_data_matrix_2, final_data_matrix_2$Sample == pt)
  
  n_clones <- as.numeric(n_distinct(sub_dat$clonality_bin, na.rm = T))
  
  if(n_clones == 0){
    n_clones=NA
  }
  
  # store the results in a dataframe
  temp <- data.frame(matrix(NA, nrow = 1, ncol = 2))
  names(temp) <- c("Sample", "n_clones")
  
  temp[1,1] <- pt
  temp[1,2] <- n_clones
  
  # Add each list in the loop to a list of lists
  z <- z + 1
  temp_list[[z]] <- temp
  
}
temp_final = as.data.frame(do.call(rbind, temp_list))
temp_final$Sample <- as.character(temp_final$Sample)
temp_final$n_clones <- as.numeric(temp_final$n_clones)

final_data_matrix_2 <- left_join(final_data_matrix_2, temp_final, by = "Sample")






# VAF correlation between mutations ####
final_data_matrix_2_sub <- subset(final_data_matrix_2, !Cohort %in% c("Au","Wang","Garg","Huet"))
final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, Gene!="MLL")

final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$mut_freq_gene > 50 & final_data_matrix_2_sub$Subset == "de_novo")

final_data_matrix_2_sub = data.frame(select(final_data_matrix_2_sub, Sample, Gene, VAF))

final_data_matrix_2_sub = final_data_matrix_2_sub %>% 
  group_by(Sample,Gene) %>% 
  summarise(VAF = max(VAF))

final_data_matrix_2_sub = dcast(final_data_matrix_2_sub, Sample~Gene, value.var="VAF")

rownames(final_data_matrix_2_sub) = final_data_matrix_2_sub$Sample
final_data_matrix_2_sub = final_data_matrix_2_sub[,-1]

final_data_matrix_2_sub[is.na(final_data_matrix_2_sub)] <- 0

corr_mat=cor(final_data_matrix_2_sub,method="p")

corrplot(corr_mat,title = "VAF Correlation", method = "square", outline = T, addgrid.col = "darkgray", order="hclust", mar = c(4,0,4,0), addrect = 8, rect.col = "black", rect.lwd = 5, cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5)





 # Age distribution by number of clones ####
 # distribution of age, etc. for clonal groups
 sub=na.omit(as.data.frame(distinct(final_data_matrix_2, Sample, Age, n_clones)))
 sub$Age=as.numeric(sub$Age)
 sub$n_clones=as.factor(sub$n_clones)


ggplot(sub, aes(x = Age, y = as.factor(n_clones), fill = n_clones)) +
  geom_density_ridges(
    jittered_points = T, position = "raincloud",
    alpha = 1, scale = 1
  ) +
   ylab(label = "Number of clones") +
   theme(legend.position = "none") +
  scale_fill_manual(name = "", values = c("1" = "#374E55FF","2" = "#DF8F44FF","3" = "#B24745FF","4" = "#79AF97FF"))

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/age_distribution_by_n_clones.png", dpi = 300, width = 5, height = 5, units = "in")
 
 
 
 
 # Mutation count by clonality ####
 # plot the correlation of mutation burden and clonal assignments
 final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
 
 final_data_matrix_2_sub <- na.omit(distinct(final_data_matrix_2_sub, Sample, n_clones, mut_freq_pt))
 final_data_matrix_2_sub$n_clones = as.factor(final_data_matrix_2_sub$n_clones)

 ggplot(final_data_matrix_2_sub, aes(x=n_clones, y=mut_freq_pt, fill = n_clones)) + 
   geom_boxplot(notch=F, outlier.colour = "white") +
   geom_jitter(shape=21, position=position_jitter(0.2)) +
   stat_compare_means() +
   theme_cowplot(font_size = 12) +
   ylim(0,11) +
   # scale_fill_uchicago() +
   scale_fill_manual(name = "", values = c("1" = "#374E55FF","2" = "#DF8F44FF","3" = "#B24745FF","4" = "#79AF97FF"))
   labs(title = NULL) +
   ylab(label= "Mutations/patient") +
   xlab(label = "Number of clones") +
     theme(legend.title = element_blank(), legend.position = "none")
   
 ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/frequency_of_mutations_in_clones_de_novo.pdf", dpi = 300, width = 4.5, height = 3, units = "in")
 
 
 # plot the distribution of number of clones per patient
 final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
 
 final_data_matrix_2_sub <- na.omit(distinct(final_data_matrix_2_sub, Sample, n_clones))
 
 ggplot(final_data_matrix_2_sub, aes(x=n_clones)) +
   geom_bar(aes(fill = as.factor(n_clones)), color = "black") +
   scale_x_discrete(limits=c("1", "2", "3", "4")) +
   scale_fill_manual(values = c("1" = "#374E55FF", "2" = "#DF8F44FF", "3" = "#B24745FF", "4" = "#79AF97FF")) +
   # scale_fill_uchicago() +
   scale_y_continuous(expand = c(0,0)) + 
   xlab("Number of clones") + 
   ylab("Number of patients") +
   theme_cowplot(font_size = 15) +
   theme(legend.title = element_blank(), legend.position = "none")
 
 
 ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/frequency_of_clones_de_novo.pdf", dpi = 300, width = 4, height = 3, units = "in")
 
 

 
 
 # survival by number of clones ####
 # final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
 final_data_matrix_2_sub = subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")

 final_data_matrix_2_sub <- final_data_matrix_2_sub[!duplicated(final_data_matrix_2_sub[1]),]
 final_data_matrix_2_sub$Time_to_OS <- (final_data_matrix_2_sub$Time_to_OS/365)
 final_data_matrix_2_sub$Censor <- as.numeric(final_data_matrix_2_sub$Censor)
 
 
 ## Kaplan Meier analysis
 # create the survival data object
 final_data_matrix_2_sub$OS <- with(final_data_matrix_2_sub, Surv(Time_to_OS, Censor ==1))
 
 # create the survival objects used to plot kaplan-meyer curves
 OS <- survfit(OS ~ n_clones, data = final_data_matrix_2_sub, conf.type = "log-log")
 
 # find the different p-values for the different comparisons
 res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ n_clones,
                          data = final_data_matrix_2_sub)
 print(res)
 
 # cohorts <- c("1", "2", "3", "4")
 # cohorts <- c("low", "high")
 
 
 # plots the survival
 surv_plot <- ggsurvplot(OS,
                         data = final_data_matrix_2_sub,
                         log = (OS),
                         log.rank.weights = c("survdiff"),
                         pval = F,
                         test.for.trend = F,
                         pval.method.size = 3,
                         pval.coord = c(0, 0),
                         conf.int = F,
                         censor = T,
                         surv.mean.line = "none",
                         risk.table = ,
                         risk.table.title = "",
                         risk.table.fontsize = 5,
                         risk.table.height = .3,
                         risk.table.y.text = T,
                         break.time.by = 5,
                         risk.table.pos = c("out"),
                         palette = c("#374E55FF", "#DF8F44FF", "#B24745FF", "#79AF97FF"),
                         xlab = "Years",
                         ylim = c(0, 1.0),
                         ylab =  "Survival Probability",
                         font.main = c(10, "plain", "#252525"),
                         pval.size = 4,
                         font.x = c(10, "plain", "#252525"),
                         font.y =  c(10, "plain", "#252525"),
                         font.legend = c(10, "plain"),
                         font.tickslab = c(10, "plain", "#252525"),
                         # legend.labs = cohorts,
                         legend.title = "Number of\nclones",
                         legend = "right",
                         ggtheme = theme(plot.title = element_text(hjust = 0.5)))
 
 print(surv_plot)
 
 png(filename = "~/Desktop/MetaAML_results/Data/Figures/meta_aml_survival_by_clone_number_de_novo.png", res = 300, width = 4, height = 2.5, units = "in")
 
 surv_plot
 print(surv_plot)
 dev.off()
 
 
 
 
 # Function: survival by mutation clonality ####
 # function to perform survival analysis for each mutation based on its clonal status
 final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
 
 final_data_matrix_2_sub$Time_to_OS <- (final_data_matrix_2_sub$Time_to_OS/365)
 final_data_matrix_2_sub$Censor <- as.numeric(final_data_matrix_2_sub$Censor)
 
 mutation_clonality_survival_function <- function(gene, range, bin, save_plot){
   
   final_data_matrix_2_sub2 <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$Gene == gene)
   final_data_matrix_2_sub2 <- final_data_matrix_2_sub2[!duplicated(final_data_matrix_2_sub2[1]),]
   
   # create the survival data object
   final_data_matrix_2_sub2$OS <- with(final_data_matrix_2_sub2, Surv(Time_to_OS, Censor ==1))
   
   if(range == T){
     # create the survival objects used to plot kaplan-meyer curves
     # final_data_matrix_2_sub2$OS <- with(final_data_matrix_2_sub2, Surv(Time_to_OS, Censor ==1))
     OS <- survfit(OS ~ clonality_bin, data = final_data_matrix_2_sub2, conf.type = "log-log")
     # find the different p-values for the different comparisons
     res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ clonality_bin,
                              data = final_data_matrix_2_sub2)
     print(res)
     cohorts <- c("< 0.2", "0.2-0.4", "0.4-0.6", "> 0.6")
     title <- "VAF group"
   }
   if(bin == T){
     # create the survival objects used to plot kaplan-meyer curves
     # final_data_matrix_2_sub2$OS <- with(final_data_matrix_2_sub2, Surv(Time_to_OS, Censor ==1))
     
     OS <- survfit(OS ~ clonality, data = final_data_matrix_2_sub2, conf.type = "log-log")
     # find the different p-values for the different comparisons
     res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ clonality,
                              data = final_data_matrix_2_sub2)
     print(res)
     # cohorts <- c("Clonal", "Clonal LOH", "Subclonal", "Subclonal LOH")
     cohorts <- c("Major", "Minor")
     
     title <- "Clonality bin"
   }
   
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
                           surv.mean.line = "none",
                           risk.table = F,
                           risk.table.title = "",
                           risk.table.fontsize = 5,
                           risk.table.height = .3,
                           risk.table.y.text = T,
                           break.time.by = 5,
                           risk.table.pos = c("out"),
                           palette = c("< 0.2" = "#374E55FF", "0.2-0.4" = "#DF8F44FF", "0.4-0.6" = "#B24745FF",  "> 0.6" = "#79AF97FF", "Major" = "#BC3C29FF", "Clonal LOH"=  "#0072B5FF", "Minor" = "#E18727FF", "Subclonal LOH" = "#20854EFF"),
                           xlab = "Years",
                           ylim = c(0, 1.0),
                           ylab =  "Survival Probability",
                           font.main = c(10, "plain", "#252525"),
                           pval.size = 4,
                           font.x = c(10, "plain", "#252525"),
                           font.y =  c(10, "plain", "#252525"),
                           font.legend = c(10, "plain"),
                           font.tickslab = c(10, "plain", "#252525"),
                           legend.labs = cohorts,
                           legend.title = paste(title),
                           legend = "right",
                           title = gene,
                           ggtheme = theme(plot.title = element_text(hjust = 0.5)))
   
   print(surv_plot)
   
   if(save_plot == T){
     png(filename = paste("~/Desktop/MetaAML_results/Data/Figures/",gene,"survival_by_VAF.png"), res = 300, width = 15, height = 5, units = "in")
     
     surv_plot
     print(surv_plot)
     dev.off()
   }
 }
 
 mutation_clonality_survival_function("NRAS", range = T, bin = F, save_plot =F)
 
 
 
 
 
# # cluster based on the VAFs ####
# final_data_matrix_2_sub <- final_data_matrix_2
# 
# final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, mut_freq_gene > 50 & Subset == "de_novo")
# final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, Gene != "MLL")
# 
# final_data_matrix_2_sub <- select(final_data_matrix_2_sub, Sample, Gene, VAF)
# final_data_matrix_2_sub$VAF <- as.numeric(as.character(final_data_matrix_2_sub$VAF))
# final_data_matrix_2_sub$VAF <- round(final_data_matrix_2_sub$VAF, 3)
# # dup <- as.data.frame(duplicated(final_data_matrix_2_sub[c("Sample", "Gene")]))
# 
# # identify patients where the same gene is mutated twice and generate a duplicate patient to represent this additional mutation
# for(i in 1:nrow(final_data_matrix_2_sub)){
#   n <- as.character(final_data_matrix_2_sub$Sample[i])
#   m <- as.character(final_data_matrix_2_sub$Gene[i])
#   sub1 <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$Sample == n & final_data_matrix_2_sub$Gene == m)
#   if(nrow(sub1) > 1){
#     for(j in 2:nrow(sub1)){
#       rn <- row_number(sub1)
#       final_data_matrix_2_sub$Sample[i] <- paste(n, "_", rn, sep="")
#     }
#   }
# }
# 
# final_data_matrix_2_sub <- dcast(final_data_matrix_2_sub, Sample ~ Gene, value.var="VAF")
# 
# rownames(final_data_matrix_2_sub) <- final_data_matrix_2_sub$Sample
# final_data_matrix_2_sub$Sample <- NULL
# final_data_matrix_2_sub[is.na(final_data_matrix_2_sub)] <- 0
# 
# # cluster and viauslize the data
# library(dendextend)
# row_dend = as.dendrogram(hclust(method = "complete", dist(final_data_matrix_2_sub)))
# row_dend = color_branches(row_dend, k = 10)
# 
# # define colors for the heatmap
# colors<-colorRampPalette(c("white","darkblue"))
# 
# 
# png(file = "~/Desktop/MetaAML_results/Data/Figures/cluster_by_vaf_2.png", res = 300, width = 5, height = 5, units = "in")
# 
# Heatmap(final_data_matrix_2_sub,
#         col = colors(10),
#         show_row_names=FALSE,
#         name = "VAF",
#         cluster_rows = T,
#         cluster_columns = T)
# dev.off()
# 
# 
# 
# library(corrplot)
# t_pts <- as.data.frame(t(final_data_matrix_2_sub))
# temp_pts <- cor(t_pts)
# 
# col<- colorRampPalette(c("darkblue", "white", "darkred"))(50)
# 
# # correlate by patients
# png(file = "~/Desktop/MetaAML_results/Data/Figures/similarity_cluster_by_patient_vaf.png", res = 300, width = 5, height = 5, units = "in")
# heatmap(x = temp_pts, col = col, symm = TRUE, labCol =FALSE, labRow =FALSE)
# dev.off()







# mutation aquizition ordering ####
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 50 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")

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

# remove FLT3 cases without VAF or defined ITD or TKD
sub=subset(sub, sub$Gene != "FLT3")

genes <- unique(sub$Gene)

# calculate frequency of gene 1 occuring before gene 2
temp_dat <- data.frame(matrix(NA, nrow = length(genes), ncol = length(genes)))

rownames(temp_dat) <- genes
colnames(temp_dat) <- genes

for(i in 1:nrow(temp_dat)){
  # select mutations of interest
  gene_x <- rownames(temp_dat)[i]
  
  for(j in 1:ncol(temp_dat)){
    gene_y <- colnames(temp_dat)[j]
    
    if(gene_x != gene_y){
      sub_gene_1 <- subset(sub, sub$Gene == gene_x)
      sub_gene_2 <- subset(sub, sub$Gene == gene_y)
      
      # select patients with both mutations
      gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
      
      n_cases <- as.numeric(nrow(gene_1_and_2))
      
      # add point color columns for visualizing clonal/subclonal trends
      gene_1_and_2$vaf_ratio <- as.numeric((gene_1_and_2$VAF_male_x.x - gene_1_and_2$VAF_male_x.y))
      
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
        
        t <- (n_1_before_2+n_2_before_1)
        
        temp_dat[i,j] <- (n_1_before_2/t) 
      }
    }
  }
}


# Get lower triangle of the correlation matrix
get_upper_tri<-function(temp_dat){
  temp_dat[lower.tri(temp_dat)] <- NA
  return(temp_dat)
}

temp_dat_final <- get_upper_tri(temp_dat)

temp_dat_final$genes <- rownames(temp_dat_final)

temp_dat_final_melted_1 <- melt(temp_dat_final, id.vars = "genes", na.rm = TRUE)

colnames(temp_dat_final_melted_1) <- c("Gene_1", "Gene_2", "fraction_1_then_2")

get_lower_tri<-function(temp_dat){
  temp_dat[upper.tri(temp_dat)] <- NA
  return(temp_dat)
}

temp_dat_final <- get_lower_tri(temp_dat)

temp_dat_final$genes <- rownames(temp_dat_final)

temp_dat_final_melted_2 <- melt(temp_dat_final, id.vars = "genes", na.rm = TRUE)

colnames(temp_dat_final_melted_2) <- c("Gene_1", "Gene_2", "fraction_1_then_2")


# combine the data together
temp_dat_final_melted <- rbind(temp_dat_final_melted_1, temp_dat_final_melted_2)



# calculate the number of cases for co-occuring mutations
sub <- subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 50 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")

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

genes <- unique(sub$Gene)

# calculate frequency of gene 1 occuring before gene 2
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
        gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF_male_x.x - gene_1_and_2$VAF_male_x.y)
        
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
temp_dat_final_melted <- left_join(temp_dat_final_melted, temp_dat_final_melted_2, by = c("Gene_1", "Gene_2"))

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

a=ggplot(data = temp_dat_final_melted, aes(y=forcats::fct_rev(reorder(Gene_1,Gene_1)), x = Gene_2)) +
  geom_point(aes(size = temp_dat_final_melted$number_1_and_2.x, color = fraction_1_then_2), shape = 15, stat = "identity") +
  geom_point(aes(size = temp_dat_final_melted$number_1_and_2.x), shape = 0, color = "#374E55FF") +
  scale_color_gradient2(high = "#003c30", low = "#543005", mid = "#f5f5f5",
                       midpoint = 0.5,
                       name= "Fraction of \nunambiguous\ncases\n(Gene 1 > Gene2)") +
  scale_size(range = c(1,7), name= "Number of\nco-occurences", breaks = c(25, 50,100,200)) +
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

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/order_of_mutations_de_novo.pdf", width = 7.5, height = 6, units = "in")






# mutation order by president vaf plot ####
# calculate the mean fraction and confidence interval for each mutation in the pairwise presidence plot by the Bradly-Terry method
# install.packages("BradleyTerryScalable")
library(BradleyTerryScalable)
# install.packages("Matrix.utils")
library(Matrix.utils)


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
      gene_1_and_2$vaf_ratio <- as.numeric((gene_1_and_2$VAF_male_x.x - gene_1_and_2$VAF_male_x.y))
      
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

# ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/bradley_terry_order_de_novo.pdf", width = 7.5, height = 6, units = "in")








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

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/vaf_distribution_order_by_time.png", dpi = 300, width = 5, height = 5, units = "in")

 

ggarrange(a,c,b,
          ncol = 3, nrow = 1, widths = c(4,1.5, 2))

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/pairwise_ordering.png", dpi = 300, width = 17, height = 7.5, units = "in")



  
  
  
  
  
  # Cox regression analysis by clonality ####
  # forrest plot
  final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo" & final_data_matrix_2$Risk != "Unknown" & final_data_matrix_2$Cohort != "Huet")
  
  # create a discrete age bin at 60 years
  final_data_matrix_2_sub$age_bin <- ifelse(final_data_matrix_2_sub$Age >= 60, 1,0)
  
  final_data_matrix_2_sub$Time_to_OS <- (final_data_matrix_2_sub$Time_to_OS/365)
  
  # forrest_plot=droplevels(final_data_matrix_2_sub)
  
  forrest_plot = droplevels(as.data.frame(unique(select(final_data_matrix_2_sub, Sample, age_bin, Sex, Risk, Cohort, n_clones.x, mut_freq_bin, mut_freq_pt, tp53_status, Shannon_diversity, age_bin, Time_to_OS, Censor))))
  
  # forrest_plot=subset(forrest_plot, forrest_plot$Gene == "NF1")
  
  forrest_plot <- within(forrest_plot, {
    n_clones <- factor(n_clones.x, labels = c("1", "2", "3", "4"))
    age_bin <- factor(age_bin, labels = c("0", "1"))
    mut_freq_bin <- factor(mut_freq_bin, labels = c("1", "2", "3", "4", "5", "6", "7", "8+"))
  })
  
  forrest_plot$Risk <- factor(forrest_plot$Risk, levels = c("Favorable","Intermediate", "Adverse"))
  forrest_plot$Risk = relevel(forrest_plot$Risk, ref = "Favorable")
  forrest_plot$Cohort = relevel(forrest_plot$Cohort, ref = "TCGA")
  forrest_plot$tp53_status <- factor(forrest_plot$tp53_status, levels = c("WT","Mut"))
  forrest_plot$tp53_status = relevel(forrest_plot$tp53_status, ref = "WT")
  forrest_plot$age_bin = relevel(forrest_plot$age_bin, ref = "0")
  
  model <- coxph( Surv(Time_to_OS, Censor) ~ Risk + age_bin + n_clones,
                  data = forrest_plot )
  ggforest(model)   # plot standard forrest plot
  
  ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/forest_plot_de_novo_nclones_risk_tp53.png", dpi = 300, width = 6, height = 7.5, units = "in")
  
  
  
  
# put the results into a dataframe for customized plotting
  forest_plot_data=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                     factor_id_sep = ":", sort_by = NULL)

  
  forest_plot_data[nrow(forest_plot_data) + 4,] = list(NA)
  forest_plot_data$factor.value[1:12]=c("Shannon Diversity", "Intermediate", "Adverse", ">60", "TP53 mut", "Majeti", "Papaemmanuil", "Tyner", "Welch", "TCGA", "TP53 wt", "Favorable")
  forest_plot_data$HR[10:12]=1 
  forest_plot_data[10:12,5:6]=1

  
  forest_plot_data$factor.value <- factor(forest_plot_data$factor.value,levels = c("Majeti", "Papaemmanuil", "Tyner", "Welch", "TCGA", "TP53 mut", "TP53 wt", "Adverse", "Intermediate", "Favorable", ">60", "Shannon Diversity"))
  
  # make custom forrest plot
  ggplot(forest_plot_data) +
    geom_hline(yintercept=1, linetype="dashed", color = "black") +
    geom_pointrange(size = .75, stat = "identity", fill = "white",
                    aes(x = factor.value, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = factor.value)) +
    scale_color_manual(values = c("Shannon Diversity" = "darkred", "Majeti"= "#636363", "Papaemmanuil"= "#636363", "Tyner"= "#636363", "Welch"= "#636363", "TCGA"= "#636363", ">60" = "#636363", "Favorable" = "#636363", "Adverse"= "#636363", "Intermediate"= "#636363", "TP53 wt"="#636363", "TP53 mut"= "#636363", "1 clone"="#374E55FF", "2 clones"= "#DF8F44FF", "3 clones"= "#B24745FF", "4 clones" ="#79AF97FF"), name = "strata") +
    ylab("Hazard Ratio")+
    theme(legend.position = "none",
      axis.title.y=element_blank()) +
    coord_flip()
  
  ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/custom_forest_plot_de_novo_shannon_diversity_risk_tp53_age.png", dpi = 300, width = 5, height = 4, units = "in")
  

  
  
  
  
  
  
  
  # vaf forrest plot ####
  # forrest plot for HRs for vafs of each mutation
  final_data_matrix_2_sub = subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 50 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")
  
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
    mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
    mut_pts = distinct(mut_pts, Sample, Gene, VAF, VAF_male_x, Time_to_OS, Censor)
    mut_pts$hr_stratifier_vaf = 0
    
    mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_male_x),]
    mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
    
    for(j in 1:nrow(mut_pts)){
      if(!is.na(mut_pts$VAF_male_x[j])){
        if(mut_pts$VAF_male_x[j] >= 0.3){
          mut_pts$hr_stratifier_vaf[j] = 1
        } 
      }
    }
    
    # run the Cox model
    model <- coxph( Surv(Time_to_OS, Censor) ~ hr_stratifier_vaf,
                    data = mut_pts )

    forest_plot_data=cox_as_data_frame(coxphsummary = model, unmangle_dict = NULL,
                                       factor_id_sep = ":", sort_by = NULL)
    forest_plot_data$gene = genes[i,1]
    
    results_list[[n]] <- forest_plot_data
    n=n+1   
  }
  
  temp_final = as.data.frame(do.call(rbind, results_list))
  
  # temp_final = subset(temp_final, temp_final$gene != "U2AF1")

  temp_final$gene <- factor(temp_final$gene, levels = temp_final$gene[order(temp_final$HR)])
  
  temp_final$sig_color = 0
  
  for(i in 1:nrow(temp_final)){
    if(temp_final$p[i] < 0.05){
      temp_final$sig_color[i] =1
    }
  }
  
  temp_final$sig_color = as.factor(temp_final$sig_color)
  
  temp_final$fdr = p.adjust(temp_final$p, method = "fdr")

  temp_final$sig_color = as.factor(temp_final$sig_color)
  # temp_final_hr_order_15$categories = as.factor(temp_final_hr_order_15$categories)
  # temp_final_hr_order_15$HR = as.numeric(temp_final_hr_order_15$HR)
  
  temp_final = subset(temp_final, temp_final$Upper_CI != "Inf")
  
  temp_final$categories <- reorder(temp_final$categories, temp_final$HR)
  
  temp_final$p_text = NA
  for(i in 1:nrow(temp_final)){
    if(temp_final$p[i] < 0.05){
      temp_final$p_text[i] = temp_final$p[i]
    }
  }
  temp_final$q_text = NA
  for(i in 1:nrow(temp_final)){
    if(temp_final$fdr[i] < 0.15){
      temp_final$q_text[i] = temp_final$fdr[i]
    }
  }
  temp_final$p_text = round(temp_final$p_text, 3)
  temp_final$q_text = round(temp_final$q_text, 2)
  
  temp_final$p_q_text = paste("p =", temp_final$p_text, "; q =", temp_final$q_text)
  temp_final$p_text = paste("p =", temp_final$p_text)
  
  for(i in 1:nrow(temp_final)){
    if(temp_final$p[i] > 0.05){
      temp_final$p_q_text[i] = ""
    }
    if(temp_final$p[i] > 0.05){
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
  
  
  ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/gene_vaf_discrete_hr_forest_plot_de_novo_30.pdf", dpi = 300, width = 5, height = 7.5, units = "in")
  
  
  
  
  

  