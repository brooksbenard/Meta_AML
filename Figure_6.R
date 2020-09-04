# ========================================================================================================================================== #
# Figure_2.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 09/01/2020
# Description: This script analyses if the VAF of a specific mutation correlates with sensitivity or resistance to different drugs in the Beat AML study
# This script will download, process, and generate the restults as seen in Figure 6 of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# =================== #
# Load libraries ####
# =================== #
if (!require('xlsx')) install.packages('xlsx'); library('xlsx')
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
if (!require('patchwork')) install.packages('patchwork'); library('patchwork')
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('ggridges')) install.packages('ggridges'); library('ggridges')

# download data from BeatAML website
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/Majeti_Lab/MetaAML_results/41586_2018_623_MOESM3_ESM.xlsx")

# data ####
# read in drug sensitivity data
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# read in BeatAML variants for analysis data
BeatAML_variants <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7)

# read in drug response data
drug_data <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 10)

# they have one drug that is not formated like all the others
drug_data[] <- lapply(drug_data, gsub, pattern = "17-AAG (Tanespimycin)", replacement = "Tanespimycin", fixed = TRUE)


# distribution of sensitivity per drug and cohort
# subset to patients with required data types
pt_subset_1 <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$exomeSeq == "y" & BeatAML_sample_data_types$totalDrug == "y")

sample_types = select(pt_subset_1, LabId, PatientId, isRelapse, isDenovo, isTransformed)
colnames(sample_types)[1] = "lab_id"

dat = na.omit(full_join(drug_data, sample_types, by = "lab_id"))

dat$subset = "Other"
for(i in 1:nrow(dat)){
  if(dat$isDenovo[i] == "TRUE" & dat$isRelapse[i] == "FALSE"){
    dat$subset[i] = "De novo"
  }
  if(dat$isDenovo[i] == "FALSE" & dat$isRelapse[i] == "TRUE"){
    dat$subset[i] = "Relapse"
  } 
  if(dat$isTransformed[i] == "TRUE" & dat$isRelapse[i] == "FALSE"){
    dat$subset[i] = "Secondary"
  }
  if(dat$isTransformed[i] == "TRUE" & dat$isRelapse[i] == "TRUE"){
    dat$subset[i] = "Relapse"
  }
}

dat$auc = as.numeric(dat$auc)

cohort_colors = list("De novo" = "#C16622FF", 
                     "Secondary" = "#767676FF", 
                     "Relapse" = "#800000FF", 
                     "Other" = "#FFA319FF")

ggplot(dat, aes(x = auc, y = as.factor(subset), fill = subset)) +
  geom_density_ridges(
    jittered_points = F, position = "raincloud",
    alpha = 1, scale = 1,
    # quantile_lines = TRUE, quantiles = 2
  ) +
  geom_vline(xintercept=70, linetype="dashed", color = "black") +
  scale_y_discrete(expand = c(0,0)) +
  ylab(label = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values = cohort_colors)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/Supplimental/beat_aml_auc_per_subset.pdf", dpi = 300, width = 5, height = 4, units = "in")

ggplot(dat, aes(x = auc, y = reorder(as.factor(inhibitor), auc))) +
  geom_density_ridges(
    alpha = 1, scale = 1
  ) +
  geom_vline(xintercept=70, linetype="dashed", color = "black") +
  scale_y_discrete(expand = c(0,0)) +
  ylab(label = NULL) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 7.5)) +
  scale_fill_discrete("grey")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/Supplimental/beat_aml_auc_per_inhibitor.pdf", dpi = 300, width = 7.5, height = 20, units = "in")



# vaf analysis ####

## Analyse all mutation relationships to a drug by performing logistic regression between  VAF vs AUC 
# select the de novo cohort for primary analysis
pt_subset_2 <- subset(pt_subset_1, pt_subset_1$isDenovo == "TRUE" & pt_subset_1$isRelapse == "FALSE")

BeatAML_variants <- BeatAML_variants %>%
  select("labId", "symbol", "t_vaf", "short_aa_change", "chrom")

# deliniate FLT3 from FLT3-ITD
for(i in 1:nrow(BeatAML_variants)){
  # print(i)
  if(BeatAML_variants[i,2] == "FLT3" & BeatAML_variants[i,4] == "ITD"){
    BeatAML_variants[i,2] <- "FLT3-ITD"
  }
  if(BeatAML_variants[i,2] == "FLT3" & BeatAML_variants[i,4] != "ITD"){
    BeatAML_variants[i,2] <- "FLT3-TKD"
  }
}


# filter to variants present in de novo patient samples
BeatAML_variants_sub <- setDT(BeatAML_variants)[labId %chin% pt_subset_2$LabId] 
BeatAML_variants_sub <- BeatAML_variants_sub %>%
  group_by(labId, symbol) %>%
  filter(t_vaf == max(t_vaf)) %>%
  ungroup()

BeatAML_variants_sub <- BeatAML_variants_sub %>%
  distinct(labId, symbol, short_aa_change, .keep_all = TRUE)

# subset to mutations present in at least five samples
mut_table <- aggregate(data.frame(count = BeatAML_variants_sub), list(value = BeatAML_variants_sub$symbol), length)
mut_table = subset(mut_table, mut_table$count.symbol > 4)

BeatAML_variants_sub <- setDT(BeatAML_variants_sub)[symbol %chin% mut_table$value] 

## pair patient samples with drug response data
# filter to only patients with both mutation and drug screening data
colnames(BeatAML_variants_sub)[1] <- "lab_id"
drug_data <- setDT(drug_data)[lab_id %chin% BeatAML_variants_sub$lab_id] 

drug_mut <- dplyr::right_join(drug_data, BeatAML_variants_sub, by = "lab_id")

# loop through drugs and create a matrix of all drug-mutation-vaf combinations
drug_list <- list()
for(i in 1:nrow(drug_mut)){
  # print(i)
  drugmut <- as.character(drug_mut[i,1])
  drugmut <- gsub("\\s.*","",drugmut)
  drug_list[[i]] <- drugmut
}
drug_mut_list <- as.data.frame(do.call(rbind, drug_list))
names(drug_mut_list) <- "Inhibitor"

drug_mut <- cbind(drug_mut_list, drug_mut)
drug_mut$inhibitor <- NULL

drug_mut$auc <- as.numeric(drug_mut$auc)
drug_mut$t_vaf <- as.numeric(drug_mut$t_vaf)

mutations <- as.data.frame(unique(drug_mut$symbol))
inhibitors_list <- na.omit(as.data.frame(unique(drug_mut$Inhibitor)))


## Find cases where the vaf and auc correlate 
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/Resistant/")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/Sensitive/")

# list for summary plots for supplimental data
vaf_auc_sensitive_plots = list()
k = 1
vaf_auc_resistant_plots = list()
l = 1

# list for aggregated summary statistics
auc_vaf_list <- list()
z <- 1

for(i in 1:nrow(inhibitors_list)){
  print(i)
  # i <- 16
  drug_sig <- as.character(inhibitors_list[i,1])
  for(j in 1:nrow(mutations)){
    # print(j)
    # j <- 1100
    mut_sig <- as.character(mutations[j,1])
    
    drug_mut_sub <- as.data.frame(drug_mut[which(drug_mut$Inhibitor == drug_sig),])
    drug_mut_sub <- as.data.frame(drug_mut_sub[which(drug_mut_sub$symbol == mut_sig),])
    
    num_obs <- as.numeric(nrow(drug_mut_sub))
    
    
    if(num_obs >= 5){
      
      # find the delta in AUN and VAF for each case in order to filter later
      delta_vaf <- diff(range(drug_mut_sub$t_vaf))
      min_vaf <- min(drug_mut_sub$t_vaf)
      max_vaf <- max(drug_mut_sub$t_vaf)
      delta_auc <- diff(range(drug_mut_sub$auc))
      min_auc <- min(drug_mut_sub$auc)
      max_auc <- max(drug_mut_sub$auc)
      
      # store all of the relationships in a dataframe
      if(num_obs >= 5 & delta_vaf >= 0.25 & delta_auc >= 75){
        
        # fit the data to a linear regression model
        lmf <- lm(drug_mut_sub$t_vaf ~ drug_mut_sub$auc, data=drug_mut_sub)
        
        # calculate the slope of the line to determine if the relationship show resistance or sensitivity
        slp <-  as.numeric(coef(lmf)[2] )
        
        # extract the r-squared value
        lmr <- summary(lmf)$r.squared
        lmr <- round(lmr, 3)
        # extract the p-value
        lmp <- as.data.frame(summary(lmf)$coefficients[,4])
        lmp <- as.numeric(lmp[2,1])
        
        auc_vaf <- data.frame(matrix(NA, nrow = 1, ncol = 11))
        names(auc_vaf) <- c("Inhibitor", "Mutated_Gene", "Slope", "R_squared", "p_value", "VAF_range", "Min_VAF", "Max_VAF", "AUC_range", "Min_AUC", "MAX_AUC")
        
        auc_vaf[1,1] <- drug_sig
        auc_vaf[1,2] <- mut_sig
        auc_vaf[1,3] <- slp
        auc_vaf[1,4] <- lmr
        auc_vaf[1,5] <- lmp
        auc_vaf[1,6] <- delta_vaf
        auc_vaf[1,7] <- min_vaf
        auc_vaf[1,8] <- max_vaf
        auc_vaf[1,9] <- delta_auc
        auc_vaf[1,10] <- min_auc
        auc_vaf[1,11] <- max_auc
        
        # Add each list in the loop to a list of lists
        auc_vaf_list[[z]] = auc_vaf 
        
        z = z + 1
        
        # bin the plots based on their category
        if(lmr >= 0.5 && lmp <= 0.05 && slp > 0 && delta_vaf >= 0.25 && delta_auc >= 75){
          
          # scatterplots ####
          p1 = ggscatter(drug_mut_sub, x = "t_vaf", y = "auc",
                         color = "black", fill = "#4393c3", size = 5, shape = 21,
                         font.label = list(color = "black", size = 9, vjust = 0.5),
                         title = paste(drug_sig, " vs. ", mut_sig, sep = ""),
                         xlab = "VAF",
                         # xlim = c(0,1),
                         ylab = "Drug AUC") +
            # ylim = c(0,300)) +
            geom_smooth(method = "lm", color = "black", alpha = .25) +
            theme(plot.title = element_text(hjust = 0.5)) +
            stat_cor(
              aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
              label.x = 0.1,
              label.y = (max_auc + 50)
            )
          
          vaf_auc_sensitive_plots[[k]] = p1
          
          k = k + 1
          
          ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/Resistant/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
        }
        if(lmr >= 0.5 && lmp <= 0.05 && slp < 0 && delta_vaf >= 0.25 && delta_auc >= 75){
          
          # make the scatterplot
          p2 = ggscatter(drug_mut_sub, x = "t_vaf", y = "auc",
                         color = "black", fill = "#b2182b", size = 5,shape = 21,
                         font.label = list(color = "black", size = 9, vjust = 0.5),
                         title = paste(drug_sig, " vs. ", mut_sig, sep = ""),
                         xlab = "VAF",
                         # xlim = c(0,1),
                         ylab = "Drug AUC") +
            # ylim = c(0,300)) +
            geom_smooth(method = "lm", color = "black", alpha = .25) +
            theme(plot.title = element_text(hjust = 0.5)) +
            stat_cor(
              aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
              label.x = 0.1,
              label.y = (max_auc + 50)
            )
          
          vaf_auc_resistant_plots[[l]] = p2
          l = l + 1
          
          ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/Sensitive/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
        }
      }
    }
  }
}
beatAML_auc_vaf_comparison <- do.call(rbind, auc_vaf_list)
beatAML_auc_vaf_comparison$fdr_corrected <- p.adjust(beatAML_auc_vaf_comparison$p_value, method = "fdr")
beatAML_auc_vaf_comparison$z_score <- qnorm(beatAML_auc_vaf_comparison$p_value)

# write out results file
write.csv(beatAML_auc_vaf_comparison, "~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation.csv", row.names=FALSE)


# make supplimental summary plots of all significant regressions
# sensitive
sensitive = list.clean(vaf_auc_sensitive_plots)

p = cowplot::plot_grid(plotlist = sensitive, 
                       ncol = 6, align = "hv",   
                       rel_heights = c(1,1),
                       rel_widths = c(1,1))

ggsave(filename ="~/Desktop/MetaAML_results/Figure_6/Supplimental/sensitive_matrix.png", dpi = 150, width = 16, height = 8, units = "in")

# resistant
resistant = list.clean(vaf_auc_resistant_plots)

p = cowplot::plot_grid(plotlist = resistant, 
                       ncol = 6, align = "hv",   
                       rel_heights = c(1,1),
                       rel_widths = c(1,1))

ggsave(filename ="~/Desktop/MetaAML_results/Figure_6/Supplimental/resistant_matrix.png", dpi = 150, width = 18, height = 15, units = "in")


# drug families for significant trends
sig_associations = read.csv("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation.csv")
sig_associations$Inhibitor = as.character(sig_associations$Inhibitor)
sig_associations = subset(sig_associations, sig_associations$p_value < 0.05)

drug_classes = read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 11)
drug_classes[] <- lapply(drug_classes, gsub, pattern = "17-AAG (Tanespimycin)", replacement = "Tanespimycin", fixed = TRUE)
colnames(drug_classes)[1] = "Inhibitor"

drug_classes$Inhibitor = gsub("\\s.*","", drug_classes$Inhibitor)

sig_genes_drug_classes = left_join(sig_associations, drug_classes, by = "Inhibitor")




# summary plots ####
# plot heatmaps of the delta AUC and delta VAF for the significant cases
vaf_data <- subset(beatAML_auc_vaf_comparison, beatAML_auc_vaf_comparison$R_squared >= 0.5 & beatAML_auc_vaf_comparison$VAF_range >= 0.25 & beatAML_auc_vaf_comparison$AUC_range >= 75)

vaf_data$sensitivity <- NA

for(i in 1:nrow(vaf_data)){
  if(vaf_data$Slope[i] < 0){
    vaf_data$AUC_range[i] <- as.numeric(as.character(vaf_data$AUC_range[i]*(-1)))
    vaf_data$sensitivity[i] <- "Sensitive"
  } else {
    vaf_data$sensitivity[i] <- "Resistant"
  }
}

occurances <- aggregate(data.frame(count = vaf_data), list(value = vaf_data$Mutated_Gene), length)
occurances <- occurances[,1:2]
colnames(occurances)[1] <- "Mutated_Gene"

vaf_data <- dplyr::left_join(vaf_data, occurances, by = "Mutated_Gene")

# add significance column for plotting
vaf_data$star = ""
for(i in 1:nrow(vaf_data)){
  if(vaf_data$fdr_corrected[i] < 0.3){
    vaf_data$star[i] = "*"
  }
}

# plot the heatmap
p = ggplot(vaf_data, aes(reorder(Mutated_Gene, -count.Inhibitor), reorder(Inhibitor, count.Inhibitor), colour = AUC_range, size = VAF_range, label = star)) +
  geom_point() +
  scale_color_gradient2(low = "#b2182b", high = "#2166ac", mid = "#f7f7f7",
                        name=expression(Delta~"AUC")) +
  geom_point(shape = 21, color = "black") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "#f0f0f0"))+
  xlab(label= NULL) +
  ylab(label="Inhibitor") +
  labs(title = NULL) +
  geom_text(size = 7.5, color = "#525252",  hjust = 0, nudge_x = 0.25)

g = guide_legend(override.aes=list(colour="lightgrey"), expression(Delta~"VAF"))
p + guides(size = g)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/Supplimentary/drug_vaf_correlation_all_summary_plot.pdf", dpi = 300, width = 6.5, height =  12, units = "in")





# binary analysis ####
## loop through all gene/drug interactions and calculate the p-value for a t.test between mutated and non-mutated samples for each drug

dir.create("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/binary/sensitive")
dir.create("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/binary/resistant")

BeatAML_variants <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7)

BeatAML_variants_sub <- BeatAML_variants %>%
  select("labId", "symbol", "t_vaf", "short_aa_change")

# deliniate FLT3 from FLT3-ITD
for(i in 1:nrow(BeatAML_variants_sub)){
  # print(i)
  if(BeatAML_variants_sub[i,2] == "FLT3" & BeatAML_variants_sub[i,4] == "ITD"){
    BeatAML_variants_sub[i,2] <- "FLT3-ITD"
  }
  if(BeatAML_variants_sub[i,2] == "FLT3" & BeatAML_variants_sub[i,4] != "ITD"){
    BeatAML_variants_sub[i,2] <- "FLT3-TKD"
  }
}

# filter to variants present in the specified subset patient group and ensure no duplicate pts are represented
BeatAML_variants_subset <- setDT(BeatAML_variants_sub)[labId %chin% pt_subset_2$LabId] 
BeatAML_variants_subset <- BeatAML_variants_subset %>%
  group_by(labId, symbol) %>%
  filter(t_vaf == max(t_vaf)) %>%
  ungroup()

BeatAML_variants_subset <- BeatAML_variants_subset %>%
  distinct(labId, symbol, .keep_all = TRUE)


drug_mut <- as.data.frame(na.omit(drug_mut))

# create a dataframe of all the unique drugs used in the screen
drug_list <- as.data.frame(unique(drug_mut$Inhibitor))
drug_list <- na.omit(drug_list)
# 122 unique drugs

datalist_beatAML <- list()
z <- 1
for(i in 1:nrow(drug_list)){
  print(i)
  
  drug <- as.character(drug_list[i,1])
  subdat <- as.data.frame(subset(drug_mut, drug_mut$Inhibitor == drug))
  
  # filter to unique patients
  subdat_uniq <- subset(subdat, !duplicated(lab_id))
  
  for(j in 1:nrow(mut_table)){
    print(j)
    # j <- 27
    gene_name <- as.character(mut_table[j,1])
    
    # create column to deliniate if the patient has a somatic mutation of interest
    subdat_uniq2 <- subdat_uniq
    
    # find pts with a mutation in the gene of interest
    pt_mut <- BeatAML_variants_subset
    pt_mut <- setDT(pt_mut)[symbol %chin% mut_table[j,1]]
    names(pt_mut) <- c("lab_id", "mutated", "t_vaf", "short_aa_change")
    
    subdat_uniq2 <- dplyr::left_join(subdat_uniq2, pt_mut, by = "lab_id")
    subdat_uniq2[is.na(subdat_uniq2)] <- "WT"
    
    mut <- subdat_uniq2[which(subdat_uniq2$mutated == gene_name),]
    nummut <- nrow(mut)
    mut <-  as.numeric(mean(mut$auc))
    
    nmut <- subdat_uniq2[which(subdat_uniq2$mutated == "WT"),]
    nmut <-  as.numeric(mean(nmut$auc))
    
    delta <- as.numeric(mut - nmut)
    
    if(nummut >= 3){
      
      if(nummut == 3){
        p_value <- as.numeric(t.test(subdat_uniq2$auc~subdat_uniq2$mutated)$p.value)
      } 
      if(nmut > 3){
        complete_dat <- nrow(na.omit(subset(subdat_uniq2, subdat_uniq2$mutated == gene_name)))
        
        if(complete_dat > 3){
          # Shapiro-Wilk normality test
          p_mut <- with(subdat_uniq2, shapiro.test(auc[mutated == gene_name]))$p.value
          p_wt <- with(subdat_uniq2, shapiro.test(auc[mutated == "WT"]))$p.value
          
          # test variance in data using F-test
          p_var <- var.test(auc ~ mutated, data = subdat_uniq2)$p.value
          
          if(p_mut > 0.05 & p_wt > 0.05 & p_var > 0.05){
            p_value <- as.numeric(t.test(subdat_uniq2$auc~subdat_uniq2$mutated)$p.value)
          } else {
            p_value <- as.numeric(wilcox.test(subdat_uniq2$auc~subdat_uniq2$mutated)$p.value)
          }  
        }
      }
      
      # plot boxplot for distribution betweeen mut and wt samples for the drug in the loop
      
      if(p_value < 0.05){
        subdat_uniq2$mutated = ifelse(subdat_uniq2$mutated != "WT", "Mut", "WT")
        
        if(delta > 0){
          
          subdat_uniq2$mutated <- factor(subdat_uniq2$mutated , levels=c("WT", "Mut"))
          
          ggplot(subdat_uniq2, aes(x=mutated, y=auc)) + 
            geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
            geom_jitter(aes(fill = mutated), color = "black", shape = 21, position=position_jitter(0.2), size = 2) +
            scale_fill_manual(values = c("WT"="lightgrey", "Mut" = "#2166ac")) +
            stat_compare_means() +
            theme_cowplot(font_size = 10) +
            labs(title = paste(drug, " vs. ", gene_name, sep = "")) +
            ylab(label= "AUC") +
            xlab(label = NULL) +
            ylim((min(subdat_uniq2$auc)), (max(subdat_uniq2$auc) + 15)) +
            theme(legend.position="none")  
          
          ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/binary/sensitive/",drug,"_",gene_name,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
        }
        if(delta < 0){
          
          subdat_uniq2$mutated <- factor(subdat_uniq2$mutated , levels=c("WT", "Mut"))
          
          ggplot(subdat_uniq2, aes(x=mutated, y=auc)) + 
            geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
            geom_jitter(aes(fill = mutated), color = "black", shape = 21, position=position_jitter(0.2), size = 2) +
            scale_fill_manual(values = c("WT"="lightgrey", "Mut" = "#b2182b")) +
            stat_compare_means() +
            theme_cowplot(font_size = 10) +
            labs(title = paste(drug, " vs. ", gene_name, sep = "")) +
            ylab(label= "AUC") +
            xlab(label = NULL) +
            ylim((min(subdat_uniq2$auc)), (max(subdat_uniq2$auc) + 15)) +
            theme(legend.position="none")   
          
          ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/binary/resistant/",drug,"_",gene_name,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
        }
      }
      
      
      beatAML_list <- data.frame(matrix(NA, nrow = 1, ncol = 6))
      names(beatAML_list) <- c("Inhibitor", "Mutated_Gene", "mean_AUC_mut", "mean_AUC_not_mut", "auc_diff", "AUC_t_test_p_value")
      
      beatAML_list[1,1] <- drug
      beatAML_list[1,2] <- gene_name
      beatAML_list[1,3] <- mut
      beatAML_list[1,4] <- nmut
      beatAML_list[1,5] <- delta
      beatAML_list[1,6] <- p_value
      
      beatAML_list <- data.frame(beatAML_list)
      z <- z + 1
      datalist_beatAML[[z]] <- beatAML_list
    }
  }
}

# concatonate all rows into one data frame
beatAML_sensitivity_difference = do.call(rbind, datalist_beatAML)

# remove cases where there wasn't a drug-gene pair
beatAML_sensitivity_difference <- beatAML_sensitivity_difference[!grepl("NaN", beatAML_sensitivity_difference$mean_AUC_mut),]  

# perform fdr correction on p-values
beatAML_sensitivity_difference$fdr_adjusted <- p.adjust(beatAML_sensitivity_difference$AUC_t_test_p_value, method = "fdr")

# add a column deliniating point color
beatAML_sensitivity_difference$color <- NA
for(i in 1:nrow(beatAML_sensitivity_difference)){
  if(beatAML_sensitivity_difference$auc_diff[i] < 0 && beatAML_sensitivity_difference$fdr_adjusted[i] <= 0.05){
    beatAML_sensitivity_difference$color[i] <- "Sensitive"
  }
  if(beatAML_sensitivity_difference$auc_diff[i] > 0 && beatAML_sensitivity_difference$fdr_adjusted[i] <= 0.05){
    beatAML_sensitivity_difference$color[i] <- "Resistant"
  }
  if(beatAML_sensitivity_difference$fdr_adjusted[i] >= 0.05){
    beatAML_sensitivity_difference$color[i] <- "NS."
  }
}

# add annotation columns
beatAML_sensitivity_difference$point_label = NA

for(i in 1:nrow(beatAML_sensitivity_difference)){
  if(-log10(beatAML_sensitivity_difference$fdr_adjusted[i]) > 2.5 & beatAML_sensitivity_difference$color[i] == "Resistant"){
    beatAML_sensitivity_difference$point_label[i] <- with(beatAML_sensitivity_difference, paste0(beatAML_sensitivity_difference$Inhibitor[i], " + ", beatAML_sensitivity_difference$Mutated_Gene[i]))
  }
  if(-log10(beatAML_sensitivity_difference$fdr_adjusted[i]) > 3){
    beatAML_sensitivity_difference$point_label[i] <- with(beatAML_sensitivity_difference, paste0(beatAML_sensitivity_difference$Inhibitor[i], " + ", beatAML_sensitivity_difference$Mutated_Gene[i]))
  }
}

beatAML_sensitivity_difference$point_label[ beatAML_sensitivity_difference$point_label == "NA" ] <- ""



# add the number of patients with mutations in each gene
mut_table_2 <- mut_table
names(mut_table_2)[1]<-"Mutated_Gene"

beatAML_sensitivity_difference <- dplyr::left_join(beatAML_sensitivity_difference, mut_table_2, by = "Mutated_Gene")

write.csv(beatAML_sensitivity_difference, "~/Desktop/MetaAML_results/Figure_6/binary_results.csv", row.names=FALSE)

# beatAML_sensitivity_difference = read.csv("~/Desktop/MetaAML_results/Data/Figures/drug_vaf_correlation/binary_response.csv")

## plot distribution results as a volcano plot
p = ggplot(beatAML_sensitivity_difference, aes(x=beatAML_sensitivity_difference$auc_diff, y=-log10(beatAML_sensitivity_difference$fdr_adjusted), color = factor(color), size = beatAML_sensitivity_difference$count.symbol)) +
  geom_point(alpha = 0.75) + 
  scale_size_continuous(range = c(.5,10)) +
  geom_point(shape = 21, color = "black", alpha = 0.25) +
  geom_label_repel(aes(label=point_label),hjust=0, vjust=1.5, size = 3) +
  scale_colour_manual(values = c("Sensitive"= "#b2182b", "Resistant"="#2166ac", "NS."="lightgrey")) + 
  theme(legend.position="right") +
  geom_hline(yintercept = 1.30103,  linetype = "dashed", color = "grey") +
  ylab(label= "-log10(FDR) q-value") +
  xlab(label= expression(Delta~"AUC (wt - mut)")) +
  labs(title = NULL) +
  theme(plot.title = element_text(color="black", size=20))

g = guide_legend("Number of\nsamples")
# h = guide_legend("")
p + guides(size = g)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/binary_results.pdf", dpi = 300, width = 6, height = 4, units = "in")

# delete large dataframes
rm(drug_list)
rm(datalist_beatAML)



# overlap ####
# overlap of significant hits for both approaches using an upset plot
binary_results=read.csv("~/Desktop/MetaAML_results/Figure_6/binary_results.csv")
binary_results=subset(binary_results, binary_results$fdr_adjusted <= 0.05)
binary_sensitive=subset(binary_results, binary_results$auc_diff < 0)
binary_resistant=subset(binary_results, binary_results$auc_diff > 0)


vaf_results=read.csv("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation.csv")
vaf_results=subset(vaf_results, vaf_results$p_value <= 0.05 & vaf_results$R_squared >= 0.5 & vaf_results$VAF_range >= 0.25 & vaf_results$AUC_range >= 75)
vaf_sensitive=subset(vaf_results, vaf_results$Slope < 0)
vaf_resistant=subset(vaf_results, vaf_results$Slope > 0)

# find overlap for interactions
b_v_s_overlap=inner_join(binary_sensitive, vaf_sensitive, by = c("Inhibitor", "Mutated_Gene"))
b_v_r_overlap=inner_join(binary_resistant, vaf_resistant, by = c("Inhibitor", "Mutated_Gene"))


upset <- data.frame(matrix(0, nrow = 110, ncol = 4))
names(upset) <- c("binary_sensitive", "binary_resistant", "vaf_sensitive", "vaf_resistant")

for(i in 1:nrow(upset)){
  if(i <= 41){
    upset$binary_sensitive[i] = 1
  }
  if(i > 41 & i <= 62){
    upset$binary_resistant[i] = 1
  }
  if(i > 62 & i <= 92){
    upset$vaf_sensitive[i] = 1
  }
  if(i > 92){
    upset$vaf_resistant[i] = 1
  }
}

# plot and save
p<-upset(upset,
         nsets = 4,
         nintersects = NA,
         sets = c("binary_sensitive", "binary_resistant", "vaf_sensitive", "vaf_resistant"),
         order.by = "freq",
         main.bar.color = "#6A6599FF",
         sets.bar.color = "#79AF97FF",
         shade.color = c("#374E55FF", "#374E55FF"),
         shade.alpha = 1,
         mainbar.y.label = "Number of\ndrug-gene pairs",
         point.size = 2,
         text.scale = 1,
         # set_size.angles = 25, 
         set_size.numbers_size = T)
pdf(file = "~/Desktop/MetaAML_results/Figure_6/Supplimental/binary_vs_vaf_drug_gene_UpSet.pdf", width = 3, height = 3)
p
dev.off()



# vaf distribution ####

BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# subset to patients with required data types
pt_subset_1 <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$exomeSeq == "y" & BeatAML_sample_data_types$totalDrug == "y")

pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" &  BeatAML_sample_data_types$isDenovo == "TRUE")


# read in BeatAML variants for analysis data
BeatAML_variants <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7)

BeatAML_variants <- BeatAML_variants %>%
  select("labId", "symbol", "t_vaf", "variant_class", "chrom")


## filter to variants present in specific disease settings
BeatAML_variants <- setDT(BeatAML_variants)[labId %chin% pt_subset$LabId]
BeatAML_variants <- BeatAML_variants %>%
  distinct(labId, symbol, .keep_all = TRUE)

mut_table <- aggregate(data.frame(count = BeatAML_variants), list(value = BeatAML_variants$symbol), length)
mut_table <- select(mut_table, "value", "count.symbol")
mut_table <- subset(mut_table, mut_table$count.symbol >= 5)

Beat_AML_muts <- setDT(BeatAML_variants)[symbol %chin% mut_table$value] 

Beat_AML_muts[] <- lapply(Beat_AML_muts, gsub, pattern = "deletion", replacement = "Del", fixed = TRUE)
Beat_AML_muts[] <- lapply(Beat_AML_muts, gsub, pattern = "insertion", replacement = "Ins", fixed = TRUE)
Beat_AML_muts[] <- lapply(Beat_AML_muts, gsub, pattern = "tandem_duplication", replacement = "ITD", fixed = TRUE)

Beat_AML_muts$t_vaf=as.numeric(Beat_AML_muts$t_vaf)

Beat_AML_muts$symbol <- with(Beat_AML_muts, reorder(symbol, -t_vaf, median))

p = ggplot(Beat_AML_muts, aes(x=symbol, y=t_vaf)) + 
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
  geom_jitter(aes(fill = variant_class), color = "black", shape = 21, position=position_jitter(0.2), size = 2) +
  scale_fill_manual(values = c("#6A6599FF", "#DF8F44FF","#79AF97FF","#B24745FF")) +
  # geom_jitter(shape=21, position=position_jitter(0.2)) +
  theme_cowplot(font_size = 10) +
  labs(title = NULL) +
  ylab(label= "VAF") +
  xlab(label = NULL) +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))

ggpar(p, legend.title = "Variant class")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/vaf_distribution_de_novo.pdf", dpi = 300, width = 7, height = 4.5, units = "in")




# co-occuring mutations and FLT3 vaf regression ####
dir.create("~/Desktop/MetaAML_results/Figure_6/FLT3_specific/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/FLT3_specific/Sensitive")

# select the different flt3-mutated samples
flt3_itd_samples = subset(drug_mut, drug_mut$symbol == "FLT3-ITD")
flt3_tkd_samples = subset(drug_mut, drug_mut$symbol == "FLT3-TKD")

background_genotype = c("DNMT3A", "NMP1", "NRAS", "TET2", "RUNX1")

for(i in 1:length(background_genotype)){
  print(background_genotype[i])
  
  # subset to only include patients with co-occuring mutations in the background genotype
  drug_mut_sub = subset(drug_mut, drug_mut$symbol == background_genotype[i]) 
  
  # subset based on flt3-itd cases
  drug_mut_sub = setDT(drug_mut_sub)[lab_id %chin% flt3_itd_samples$lab_id]
  drug_mut_sub = select(drug_mut_sub, lab_id, symbol, t_vaf)
  names(drug_mut_sub) = c("lab_id", "gene2", "gene2_vaf")
  
  itd = inner_join(flt3_itd_samples, drug_mut_sub, by = "lab_id")
  
  # subset based on flt3-tkd cases
  drug_mut_sub = subset(drug_mut, drug_mut$symbol == background_genotype[i]) 
  drug_mut_sub = setDT(drug_mut_sub)[lab_id %chin% flt3_tkd_samples$lab_id]
  drug_mut_sub = select(drug_mut_sub, lab_id, symbol, t_vaf)
  names(drug_mut_sub) = c("lab_id", "gene2", "gene2_vaf")
  
  tkd = inner_join(flt3_tkd_samples, drug_mut_sub, by = "lab_id")
  
  data_list = list(itd,tkd)
  
  auc_vaf_list <- list()
  vaf_auc_sensitive_plots = list()
  
  z = 1
  k = 1
  
  for(a in 1:length(data_list)){
    
    sub_data = as.data.frame(data_list[a])
    flt3_type = sub_data$symbol[1]
    
    for(j in 1:nrow(inhibitors_list)){
      drug_sig <- as.character(inhibitors_list[j,1])
      
      sub = unique(subset(sub_data, sub_data$inhibitor == drug_sig))
      
      sub$auc = as.numeric(sub$auc)
      
      mut_sig = sub$gene2[1]
      
      if(nrow(sub) >= 5){
        # print(i)
        # find the delta in AUN and VAF for each case in order to filter later
        delta_vaf <- as.numeric(diff(range(sub$t_vaf)))
        min_vaf <- as.numeric(min(sub$t_vaf))
        max_vaf <- as.numeric(max(sub$t_vaf))
        delta_auc <- as.numeric(diff(range(sub$auc)))
        min_auc <- as.numeric(min(sub$auc))
        max_auc <- as.numeric(max(sub$auc))
        
        # store all of the relationships in a dataframe
        if(delta_vaf >= 0.25 & delta_auc >= 75){
          
          # fit the data to a linear regression model
          lmf <- lm(sub$t_vaf ~ sub$auc, data=sub)
          
          # calculate the slope of the line to determine if the relationship show resistance or sensitivity
          slp <-  as.numeric(coef(lmf)[2] )
          
          # extract the r-squared value
          lmr <- summary(lmf)$r.squared
          lmr <- round(lmr, 3)
          # extract the p-value
          lmp <- as.data.frame(summary(lmf)$coefficients[,4])
          lmp <- as.numeric(lmp[2,1])
          
          auc_vaf <- data.frame(matrix(NA, nrow = 1, ncol = 11))
          names(auc_vaf) <- c("Inhibitor", "Mutated_Gene", "Slope", "R_squared", "p_value", "VAF_range", "Min_VAF", "Max_VAF", "AUC_range", "Min_AUC", "MAX_AUC")
          
          auc_vaf[1,1] <- drug_sig
          auc_vaf[1,2] <- mut_sig 
          auc_vaf[1,3] <- slp
          auc_vaf[1,4] <- lmr
          auc_vaf[1,5] <- lmp
          auc_vaf[1,6] <- delta_vaf
          auc_vaf[1,7] <- min_vaf
          auc_vaf[1,8] <- max_vaf
          auc_vaf[1,9] <- delta_auc
          auc_vaf[1,10] <- min_auc
          auc_vaf[1,11] <- max_auc
          
          # Add each list in the loop to a list of lists
          auc_vaf_list[[z]] = auc_vaf 
          
          z = z + 1
          
          # bin the plots based on their category
          if(lmr >= 0.5 && lmp <= 0.05 && slp > 0 && delta_vaf >= 0.25 && delta_auc >= 75){
            
            # scatterplots ####
            p1 = ggscatter(sub, x = "t_vaf", y = "auc",
                           color = "black", fill = "#4393c3", size = 5, shape = 21,
                           font.label = list(color = "black", size = 9, vjust = 0.5),
                           title = paste(drug_sig, "\n", mut_sig, " + ", flt3_type, sep = ""),
                           xlab = paste(flt3_type, "VAF", sep = " "),
                           ylab = "Drug AUC") +
              geom_smooth(method = "lm", color = "black", alpha = .25) +
              theme(plot.title = element_text(hjust = 0.5)) +
              stat_cor(
                aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                label.x = 0.1,
                label.y = (max_auc + 50)
              )
            
            vaf_auc_sensitive_plots[[k]] = p1
            
            k = k + 1
            
            ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/FLT3_specific/Resistant/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
          }
          if(lmr >= 0.5 && lmp <= 0.05 && slp < 0 && delta_vaf >= 0.25 && delta_auc >= 75){
            #   
            # make the scatterplot
            p2 = ggscatter(sub, x = "t_vaf", y = "auc",
                           color = "black", fill = "#b2182b", size = 5,shape = 21,
                           font.label = list(color = "black", size = 9, vjust = 0.5),
                           title = paste(drug_sig, "\n", mut_sig, " + ", flt3_type, sep = ""),
                           xlab = paste(flt3_type, "VAF", sep = " "),
                           # xlim = c(0,1),
                           ylab = "Drug AUC") +
              # ylim = c(0,300)) +
              geom_smooth(method = "lm", color = "black", alpha = .25) +
              theme(plot.title = element_text(hjust = 0.5)) +
              stat_cor(
                aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                label.x = 0.1,
                label.y = (max_auc + 50)
              )
            
            vaf_auc_resistant_plots[[l]] = p2
            l = l + 1
            
            ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/FLT3_specific/Sensitive/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
          }
        }
      }
    }  
  }
  
}
# 
flt3_auc_vaf_comparison <- do.call(rbind, auc_vaf_list)
flt3_auc_vaf_comparison$fdr_corrected <- p.adjust(flt3_auc_vaf_comparison$p_value, method = "fdr")
flt3_auc_vaf_comparison$z_score <- qnorm(flt3_auc_vaf_comparison$p_value)

# write out results file
write.csv(flt3_auc_vaf_comparison, "~/Desktop/MetaAML_results/Figure_6/FLT3_specific/flt3_drug_vaf_correlation.csv", row.names=FALSE)
