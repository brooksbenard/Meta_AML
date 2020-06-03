#### load required packages and data ####
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
#
## source data compile script
# source("~/Desktop/MetaAML_results/MetaAML_data_compile_script.R")
# load R data file from data compile script
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# source("~/Desktop/Majeti_Lab/Scripts/TCGA_BeatAML_Stanford_mutation_aggrigation.R")
# 
# AML_pooled_co_mutation_function(pt_subset = "de novo", vaf_or_variant = "variant", n_muts = 5, subset_to_genes = F, include = F, gene_list = c("KRAS","NRAS"), file_label = "KRAS_NRAS_included", survival_analysis = F, survival_mutations = c("WT1", "NRAS"), grouped_or_individual = "individual", save_files = T)


MetaAML_survival_co_mutation_analysis_function <- function(pt_subset, n_muts){
  # subset to desired cohorts
  if(pt_subset == "De novo" | pt_subset == "De_novo"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
    label1 <- "De novo"
    label2 <- "De_novo"
    
  }
  if(pt_subset == "Secondary"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "transformed")
    label1 <- "Secondary"
    label2 <- "Transformed"
    
  }
  if(pt_subset == "Relapse"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
    label1 <- "Relapse"
    label2 <- "Relapse"
    
  }

  
  ### recalculate the frequency of mutations per gene and number of mutations per patient sample
  # find the frequency of mutations per patient
  combined_freq <- aggregate(data.frame(count = final_data_matrix_sub), list(value = final_data_matrix_sub$Sample), length)
  
  combined_pt_freq <- select(combined_freq, "value", "count.Sample")
  colnames(combined_pt_freq)[1] <- "Sample"
  
  combined_mid <- dplyr::left_join(final_data_matrix_sub, combined_pt_freq, by = "Sample")
  
  # find the frequency of mutations per gene
  combined_mut_freq <- aggregate(data.frame(count = final_data_matrix_sub), list(value = final_data_matrix_sub$symbol), length)
  
  combined_mut_freq <- select(combined_mut_freq, "value", "count.symbol")
  colnames(combined_mut_freq)[1] <- "symbol"
  
  combined_final <- dplyr::left_join(combined_mid, combined_mut_freq, by = "symbol")
  
  
  ## for aesthetic purposes, filter to only genes mutated at least "n" times
  n <- as.numeric(n_muts)
  
  combined_final_sub <- subset(combined_final, combined_final$count.symbol >= n)
  
  
  
  # make list of recurrent mutations
  mutations <- as.data.frame(unique(as.character(combined_final_sub$symbol)))
  colnames(mutations)[1] <- "symbol"
  mutations$symbol <- as.character(mutations$symbol)
  
  n_muts_1 <- as.numeric(nrow(mutations))
  n_muts_2 <- (n_muts_1 + 1)
  
  # make list of patient ids
  patients <- as.data.frame(unique(combined_final_sub$Sample))
  colnames(patients)[1] <- "Sample"
  
  n_pts <- as.numeric(nrow(patients))
  
  # make new data frame to populate with results
  temp_dat2 <- data.frame(matrix(NA, nrow = n_pts, ncol = n_muts_2))
  colnames(temp_dat2)[2:n_muts_2] <- c(mutations$symbol)
  
  colnames(temp_dat2)[1] <- "Sample"
  temp_dat2$Sample <- as.character(patients$Sample)
  
  n_pts <- as.numeric(nrow(temp_dat2))
  
  for(i in 1:nrow(combined_final_sub)){
    # i <- 2
    id <- as.character(combined_final_sub$Sample[i])
    mut <- as.character(combined_final_sub$symbol[i])
    
    r_num <- as.numeric(which(temp_dat2$Sample == id))
    c_num <- as.numeric(which(colnames(temp_dat2) == mut))
    
    temp_dat2[r_num,c_num] <- 1
  }
  
  temp_dat2[is.na(temp_dat2)] <- 0
  
  temp_dat2 <- as.data.frame(temp_dat2)
  
  combined_final_sub <- combined_final_sub[order(combined_final_sub$count.symbol, decreasing = T),]
  
  
  sort_desc <- c(as.character(combined_final_sub$symbol))
  
  
  # order the plot
  temp_dat2 <- temp_dat2 %>% 
    arrange_at(sort_desc, desc)
  
  
  pt_order <- as.data.frame(temp_dat2$Sample)
  colnames(pt_order)[1] <- "Sample"
  pt_order$order <- seq.int(nrow(pt_order))
  
  mut_table_final <- left_join(pt_order, combined_final_sub, by = "Sample")
  
  mut_table_final$variant_type <- as.character(mut_table_final$variant_type)
  
  # manually annotate 'insertion' and 'tandem duplications' variants as ITD in FLT3 cases
  for (i in 1:nrow(mut_table_final)) {
    if(mut_table_final$symbol[i] == "FLT3" & mut_table_final$variant_type[i] == "insertion"){
      mut_table_final$variant_type[i] <- "ITD"
    }
    if(mut_table_final$variant_type[i] == "tandem_duplication"){
      mut_table_final$variant_type[i] <- "ITD"
    }
  }
  
  
  #### perform Fisher's exact test on mutation occurences in the combined data set ####
  # deliniate between FLT3-ITD and FLT3-TKD because they will have different patterns with mutations (...I think)
  
  # create fresh dataframe
  final_data_matrix <- mut_table_final
  
  # manually annotate FLT as TKD or ITD
  for (i in 1:nrow(final_data_matrix)) {
    if(final_data_matrix$symbol[i] == "FLT3" & final_data_matrix$variant_type[i] == "ITD"){
      final_data_matrix$symbol[i] <- "FLT3_ITD"
    }
    if(final_data_matrix$symbol[i] == "FLT3" & final_data_matrix$variant_type[i] == "deletion"){
      final_data_matrix$symbol[i] <- "FLT3_ITD"
    }
    if(final_data_matrix$symbol[i] == "FLT3" & final_data_matrix$variant_type[i] == "SNV"){
      final_data_matrix$symbol[i] <- "FLT3"
    }
  }
  
  # make list of recurrent mutations
  mutations_2 <- as.data.frame(unique(as.character(final_data_matrix$symbol)))
  colnames(mutations_2)[1] <- "symbol"
  mutations_2$symbol <- as.character(mutations_2$symbol)
  
  n_muts_2 <- as.numeric(nrow(mutations_2))
  
  n_muts_3 <- (n_muts_2 + 1)
  
  # make new data frame to populate
  temp_dat3 <- data.frame(matrix(NA, nrow = n_pts, ncol = n_muts_3))
  colnames(temp_dat3)[2:n_muts_3] <- c(mutations_2$symbol)
  
  temp_dat3[,order(colnames(temp_dat3))]
  
  colnames(temp_dat3)[1] <- "Sample"
  temp_dat3$Sample <- as.character(patients$Sample)
  
  temp_dat3 <- temp_dat3[,order(colnames(temp_dat3))]
  
  n_pts <- as.numeric(nrow(temp_dat3))
  
  for(i in 1:nrow(final_data_matrix)){
    # i <- 2
    id <- as.character(final_data_matrix$Sample[i])
    mut <- as.character(final_data_matrix$symbol[i])
    
    r_num <- as.numeric(which(temp_dat3$Sample == id))
    c_num <- as.numeric(which(colnames(temp_dat3) == mut))
    
    temp_dat3[r_num,c_num] <- 1
  }
  
  temp_dat3[is.na(temp_dat3)] <- 0
  
  temp_dat3 <- as.data.frame(temp_dat3)
  
  
  
  # create new dataframe for test results
  f_test_matrix <- data.frame(matrix(NA, nrow = n_muts_2, ncol = n_muts_2))
  temp_dat3$Sample <- NULL
  c_names <- colnames(temp_dat3)
  
  colnames(f_test_matrix) <- c(c_names)
  rownames(f_test_matrix) <- c(c_names)
  
  # create dataframe for p-values from fisher's exact test
  f_test_matrix_p <- f_test_matrix
  
  for(i in 1:nrow(f_test_matrix)){
    # i <- 1
    r_name <- as.character(rownames(f_test_matrix[i,]))
    for(j in 1:ncol(f_test_matrix)){
      if(i < j | i > j){
        # j <- 10
        c_name <- as.character(colnames(f_test_matrix)[j])
        
        # select the columns in the binarized mutation file
        c2_num <- as.numeric(grep(r_name, colnames(temp_dat3)))
        
        if(r_name == "FLT3" | c_name == "FLT3"){
          c2_num <- c2_num[1]
        }
        
        c1_num <- as.numeric(grep(c_name, names(temp_dat3)))
        
        if(r_name == "FLT3" | c_name == "FLT3"){
          c1_num <- c1_num[1]
        }
        
        f_test <- as.data.frame(temp_dat3[,c(c1_num,c2_num)])
        f_test$sum <- (f_test[,1] + f_test[,2])
        
        f_test_sub <- data.frame(matrix(NA, nrow = 2, ncol = 2))
        
        colnames(f_test_sub)[1:2] <- c(r_name, c_name)
        rownames(f_test_sub)[1:2] <- c(c_name, r_name)
        
        f_test_sub[1,1] <- as.numeric(nrow(subset(f_test, f_test[,1] == 0 & f_test[,2] == 0)))
        f_test_sub[1,2] <- as.numeric(nrow(subset(f_test, f_test[,1] == 1)))
        f_test_sub[2,1] <- as.numeric(nrow(subset(f_test, f_test[,2] == 1)))
        f_test_sub[2,2] <- as.numeric(nrow(subset(f_test, f_test[,1] == 1 & f_test[,2] == 1)))
        
        diff <- (f_test_sub[2,2])-(((f_test_sub[1,2])/n_pts)*((f_test_sub[2,1])/n_pts)*n_pts)
        # fold_change <- (f_test_sub[2,2])/(((f_test_sub[1,2])/n_pts)*((f_test_sub[2,1])/n_pts)*n_pts)
        
        
        p_val <- fisher.test(f_test_sub)$p.value
        
        f_test_matrix[i,j] <- diff
        f_test_matrix_p[i,j] <- p_val 
      }
    }
  }
  
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(f_test_matrix){
    f_test_matrix[lower.tri(f_test_matrix)] <- NA
    return(f_test_matrix)
  }
  
  final_cor_frame_d <- get_lower_tri(f_test_matrix)
  
  final_cor_frame_d$genes <- rownames(final_cor_frame_d)
  
  # Melt the correlation matrix
  melted_cormat_d <- melt(final_cor_frame_d, na.rm = TRUE)
  colnames(melted_cormat_d)[3] <- "Delta"
  
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(f_test_matrix_p){
    f_test_matrix_p[lower.tri(f_test_matrix_p)] <- NA
    return(f_test_matrix_p)
  }
  
  final_cor_frame_p <- get_lower_tri(f_test_matrix_p)
  
  final_cor_frame_p$genes <- rownames(final_cor_frame_p)
  
  # Melt the correlation matrix
  melted_cormat_p <- melt(final_cor_frame_p, na.rm = TRUE)
  colnames(melted_cormat_p)[3] <- "P_value"
  
  # create the final dataframe with all variables used for plotting
  melted_cormat <- left_join(melted_cormat_p, melted_cormat_d, by = c("genes", "variable"))
  colnames(melted_cormat)[1:4] <- c("Gene_1", "Gene_2", "P_value", "Delta")
  
  
  # for each significantly co-occuring or mutually exclusive interaction 
  melted_cormat$P_value <- as.numeric(as.character(melted_cormat$P_value))
  
  melted_cormat$fdr_q_value <- p.adjust(melted_cormat$P_value, method = "fdr")
  
  # remove all p values that are not significant
  for(i in 1:nrow(melted_cormat)){
    if(melted_cormat$fdr_q_value[i] > 0.1){
      melted_cormat$fdr_q_value[i] <- NA
    }
  }

  
  melted_cormat_sub <- subset(melted_cormat, melted_cormat$fdr_q_value <= 0.1)
  melted_cormat_sub$Gene_2 <- as.character(melted_cormat_sub$Gene_2)
  rownames(melted_cormat_sub) <- 1:nrow(melted_cormat_sub)
  
  final_data_matrix_sub_copy <- final_data_matrix_sub
  
  for(i in 1:nrow(melted_cormat_sub)){
    print(i)
    # i <- 1
    gene1 <- melted_cormat_sub$Gene_1[i]
    gene2 <- melted_cormat_sub$Gene_2[i]
    
    # find patients with mutations in the genes
    final_data_matrix_sub1 <- subset(final_data_matrix_sub_copy, final_data_matrix_sub_copy$symbol == gene1)
    final_data_matrix_sub2 <- subset(final_data_matrix_sub_copy, final_data_matrix_sub_copy$symbol == gene2)
    final_data_matrix_sub3 <- setDT(final_data_matrix_sub1)[(Sample) %chin% final_data_matrix_sub2$Sample]
    
    if(nrow(final_data_matrix_sub1) > 1 & nrow(final_data_matrix_sub2) > 1 & nrow(final_data_matrix_sub3) > 1){
      # create strata columns
      
      final_data_matrix_sub1$STRATA1 <- 1
      final_data_matrix_sub2$STRATA2 <- 2
      final_data_matrix_sub3$STRATA3 <- 3
      
      # subset to informative columns
      final_data_matrix_sub1 <- select(final_data_matrix_sub1, Sample, STRATA1)
      final_data_matrix_sub2 <- select(final_data_matrix_sub2, Sample, STRATA2)
      final_data_matrix_sub3 <- select(final_data_matrix_sub3, Sample, STRATA3)
      
      # add strata to the larger data matrix
      final_data_matrix_sub_copy1 <- left_join(final_data_matrix_sub_copy, final_data_matrix_sub1, by = "Sample")
      final_data_matrix_sub_copy1 <- left_join(final_data_matrix_sub_copy1, final_data_matrix_sub2, by = "Sample")
      final_data_matrix_sub_copy1 <- left_join(final_data_matrix_sub_copy1, final_data_matrix_sub3, by = "Sample")
      
      # create consensus strata column
      final_data_matrix_sub_copy1$STRATA <- rowSums(final_data_matrix_sub_copy1[,c("STRATA1", "STRATA2", "STRATA3")], na.rm=TRUE)
      
      final_data_matrix_sub_copy1[,12:14] <- NULL
      
      # filter to unique patients
      final_data_matrix_sub_final <- final_data_matrix_sub_copy1[!duplicated(final_data_matrix_sub_copy1[1]),]
      final_data_matrix_sub_final <- na.omit(final_data_matrix_sub_final)
      
      nstrata <- as.numeric(length(unique(final_data_matrix_sub_final$STRATA)))
      
      if(nstrata == 4){
        c <- c("grey50", "darkblue", "darkgreen", "darkred")
        l_labs <- c("neither", gene1, gene2, "both")
        
        final_data_matrix_sub_final$Time_to_OS <- (final_data_matrix_sub_final$Time_to_OS/365)
        
        # create the survival data 
        final_data_matrix_sub_final$OS <- with(final_data_matrix_sub_final, Surv(Time_to_OS, Censor == 1))
        
        # create the survival objects used to plot kaplan-meyer curves
        OS <- survfit(OS ~ STRATA, data = final_data_matrix_sub_final, conf.type = "log-log")
        
        
        # find the different p-values for the different comparisons
        # final_data_matrix_sub_final$Censor <- as.numeric(final_data_matrix_sub_final$Censor)  
        # res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ STRATA,
        #                          data = final_data_matrix_sub_final)
        # print(res)
        
        
        # plots the survival
        surv_plot <- ggsurvplot(OS,
                                data = final_data_matrix_sub_final,
                                log = (OS),
                                log.rank.weights = c("survdiff"),
                                pval = T,
                                test.for.trend = T,
                                pval.method.size = 3,
                                pval.coord = c(0, 0),
                                conf.int = F,
                                censor = F,
                                surv.median.line = "none",
                                risk.table = T,
                                risk.table.title = "",
                                risk.table.fontsize = 5,
                                risk.table.height = .25,
                                risk.table.y.text = T,
                                break.time.by = 1,
                                risk.table.pos = c("out"),
                                palette = c,
                                title = paste("MetaAML survival by ", gene1, " and ", gene2, " mutations ", "(",label1,")", sep = ""),
                                xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(20, "plain", "black"),
                                pval.size = 5,
                                font.x = c(20, "plain", "black"),
                                font.y =  c(20, "plain", "black"),
                                font.legend = c(15, "plain"),
                                font.tickslab = c(15, "plain", "black"),
                                legend.labs = l_labs,
                                legend.title = "Mutation status",
                                legend = "right",
                                ggtheme = theme(plot.title = element_text(hjust = 0.5)))
        
        print(surv_plot)
        
        png(filename = paste("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/Survival_by_significant_cases/",label2,"/MetaAML_survival_", gene1, "_", gene2,"_", pt_subset, ".png", sep = ""), res = 300, width = 15, height = 10, units = "in")
        surv_plot
        print(surv_plot)
        dev.off()
      }
      
      
      if(nstrata == 3){
        c <- c("darkblue", "darkgreen", "darkred")
        l_labs <- c("neither", gene1, gene2)
        
        final_data_matrix_sub_final$Time_to_OS <- (final_data_matrix_sub_final$Time_to_OS/365)
        
        # create the survival data 
        final_data_matrix_sub_final$OS <- with(final_data_matrix_sub_final, Surv(Time_to_OS, Censor == 1))
        
        # create the survival objects used to plot kaplan-meyer curves
        OS <- survfit(OS ~ STRATA, data = final_data_matrix_sub_final, conf.type = "log-log")
        
        
        # find the different p-values for the different comparisons
        # final_data_matrix_sub_final$Censor <- as.numeric(final_data_matrix_sub_final$Censor)  
        # res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ STRATA,
        #                          data = final_data_matrix_sub_final)
        # print(res)
        
        
        # plots the survival
        surv_plot <- ggsurvplot(OS,
                                data = final_data_matrix_sub_final,
                                log = (OS),
                                log.rank.weights = c("survdiff"),
                                pval = T,
                                test.for.trend = T,
                                pval.method.size = 3,
                                pval.coord = c(0, 0),
                                conf.int = F,
                                censor = F,
                                surv.median.line = "none",
                                risk.table = T,
                                risk.table.title = "",
                                risk.table.fontsize = 5,
                                risk.table.height = .25,
                                risk.table.y.text = T,
                                break.time.by = 1,
                                risk.table.pos = c("out"),
                                palette = c,
                                title = paste("MetaAML survival by ", gene1, " and ", gene2, " mutations ", "(",label1,")", sep = ""),
                                xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(20, "plain", "black"),
                                pval.size = 5,
                                font.x = c(20, "plain", "black"),
                                font.y =  c(20, "plain", "black"),
                                font.legend = c(15, "plain"),
                                font.tickslab = c(15, "plain", "black"),
                                legend.labs = l_labs,
                                legend.title = "Mutation status",
                                legend = "right",
                                ggtheme = theme(plot.title = element_text(hjust = 0.5)))
        
        print(surv_plot)
        
        png(filename = paste("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/Survival_by_significant_cases/",label2,"/MetaAML_survival_", gene1, "_", gene2,"_", pt_subset, ".png", sep = ""), res = 150, width = 15, height = 10, units = "in")
        surv_plot
        print(surv_plot)
        dev.off()
      }
    }
  }
}

MetaAML_survival_co_mutation_analysis_function(pt_subset = "De_novo", n_muts = 5)
