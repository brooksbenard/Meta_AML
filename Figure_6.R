# ========================================================================================================================================== #
# Figure_6.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 11/02/2021
# Description: This script analyses if the VAF of a specific mutation correlates with sensitivity or resistance to different drugs in the Beat AML study
# This script will download, process, and generate the restults as seen in Figure 6 of the manuscript Benard et al. "Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia".
# ========================================================================================================================================== #

# ================ #
# Load packages ####
# ================ #
# Package names
packages <- c("ggplot2", "tidyverse", "viridis", "viridisLite", "RColorBrewer", "tydyr", "dplyr", "cometExactTest", "discover", "stringr", "maditr", "reshape2", "data.table", "epitools", "corrplot", "plyr", "muhaz", "reshape", "survival", "survivalAnalysis", "survMisc", "survminer", "ggsci", "vegan", "ggrepel", "ggforce", "rstatix", "effsize", "psych", "maxstat", "RCurl", "ggpubr", "UpSetR", "cowplot", "readxl", "scales", "rlist", "ggplotify", "ggridges", "gridExtra", "magrittr", "stats", "patchwork", "xlsx")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# create directores
dir.create("~/Desktop/MetaAML_results/Figure_6")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental")

# data ####
# download data from BeatAML website
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/Majeti_Lab/MetaAML_results/41586_2018_623_MOESM3_ESM.xlsx")


# VAF correction ####
# correct the VAFs based on reported cytogenetic data per patient
# read in BeatAML variants for analysis'
BeatAML_variants <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7, col_names = TRUE)

# extract useful columns
BeatAML_variants_sub <- BeatAML_variants %>%
  dplyr::select("labId", "symbol", "chrom", "t_vaf", "variant_class", "short_aa_change", "genotyper")

## Keeping only the hits that are found by both mutect and varscan. Also keeping all pindel hits
# arrange by VAF in descending order so that we keep only the highest VAF from overlapping hits
BeatAML <- BeatAML_variants_sub %>%
  arrange(desc(t_vaf))

# create new data frame with the columns needed to determine mutation overlaps 
overlaps <- BeatAML %>%
  select(labId, symbol, short_aa_change)

# create a column in BeatAML where TRUE if the value is duplicated earlier in the data table (does not put TRUE for first instance of that duplicate)
BeatAML$duplicates <- NA
BeatAML$duplicates <- duplicated(overlaps)

# create a column that lists TRUE if the row should be kept and FALSE if it should be removed
BeatAML$keep <- NA
for (i in 1:nrow(BeatAML)) {
  if (BeatAML$genotyper[i] == "pindel" | BeatAML$duplicates[i] | BeatAML$symbol[i] == "NPM1") {
    BeatAML$keep[i] <- TRUE
  } else {
    BeatAML$keep[i] <- FALSE
  }
}

# keep only those that are true in keep column
BeatAML <- BeatAML[which(BeatAML$keep),]

## change back to being arranged by LabId
BeatAML_variants_sub <- BeatAML %>%
  arrange(labId)


colnames(BeatAML_variants_sub)[1] = "LabId"

# read in cytogenetic data 
BeatAML_cytogenetics = read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# extract useful columns
BeatAML_cytogenetics_sub <- BeatAML_cytogenetics %>%
  dplyr::select("LabId", "consensus_sex", "Karyotype", "Other.Cytogenetics")

BeatAML_aggrigate = left_join(BeatAML_variants_sub, BeatAML_cytogenetics_sub, by = "LabId")

# subset to frequent mutations
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1","NRAS", "KRAS", "FLT3", "PTPN11")

BeatAML_aggrigate = subset(BeatAML_aggrigate, BeatAML_aggrigate$symbol %in% mut_list)

# download a file correlating the chromosome loci with gene ID
# http://uswest.ensembl.org/biomart/martview/a3043f537a26692111e9a7a94003ff68
all_genes <- read.table("~/Downloads/mart_export.txt", sep = "\t", header = T, fill = T,  stringsAsFactors = FALSE, quote = "")

all_genes = subset(all_genes, all_genes$Gene.name %in% mut_list)

all_genes$Karyotype_arm = gsub("\\..*","",all_genes$Karyotype.band)
all_genes$partial_annotation = paste("(",all_genes$Chromosome.scaffold.name,")","(",all_genes$Karyotype_arm,")", sep = "")
all_genes$full_annotation = paste("(",all_genes$Chromosome.scaffold.name,")","(",all_genes$Karyotype.band,")", sep = "")

# append to mutation file
names(all_genes)[1] = "symbol"

BeatAML_aggrigate = left_join(BeatAML_aggrigate, all_genes, by = "symbol")

BeatAML_aggrigate$Karyotype = gsub('\\s+', '', BeatAML_aggrigate$Karyotype)
BeatAML_aggrigate$VAF_CN_corrected = NA

BeatAML_aggrigate$t_vaf = as.numeric(BeatAML_aggrigate$t_vaf*100)

for(i in 1:nrow(BeatAML_aggrigate)){
  
  # find whole chromosome gains or losses
  ch_gain = paste("+",BeatAML_aggrigate$chrom[i], sep = "")
  ch_loss = paste("-",BeatAML_aggrigate$chrom[i], sep = "")
  
  # correct for VAF based on copy number gains for that chromosome  
  if(grepl(ch_gain, BeatAML_aggrigate$Karyotype[i], fixed = T) & BeatAML_aggrigate$t_vaf[i] > 35) {BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]*0.75}
  if(grepl(ch_gain, BeatAML_aggrigate$Karyotype[i], fixed = T) & BeatAML_aggrigate$t_vaf[i] <= 35) {BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]*1.5}
  
  # correct for VAF based on copy number loss for that chromosome
  if(grepl(ch_loss, BeatAML_aggrigate$Karyotype[i])){ BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]/2}
  
  # find focal chromosome gains or losses
  locus_gain = as.character(paste("add",BeatAML_aggrigate$partial_annotation[i],"|","add",BeatAML_aggrigate$full_annotation[i], sep = ""))
  locus_loss = as.character(paste("del",BeatAML_aggrigate$partial_annotation[i],"|","del",BeatAML_aggrigate$full_annotation[i], sep = ""))
  
  # correct for VAF based on copy number gains at the gene locus  
  if(grepl(locus_gain, BeatAML_aggrigate$Karyotype[i], fixed = T) & BeatAML_aggrigate$t_vaf[i] > 35) { BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]*0.75}
  if(grepl(locus_gain, BeatAML_aggrigate$Karyotype[i], fixed = T) & BeatAML_aggrigate$t_vaf[i] <= 35) { BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]*1.5}
  
  # correct for VAF based on copy number loss at the gene locus
  if(grepl(locus_loss, BeatAML_aggrigate$Karyotype[i], fixed = T)) { BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]/2}
  
}
# TP53 is a unique case. Code manually for focal deletions
tp53_loss1 = "del(17)(p13)"
tp53_loss2 = "del(17)(p11.2p13)"

for(i in 1:nrow(BeatAML_aggrigate)){
  # if no copy number differenes detected, populate the raw VAF
  if(grepl(tp53_loss1, BeatAML_aggrigate$Karyotype[i]) & BeatAML_aggrigate$symbol[i] == "TP53") { BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]/2}
  
  if(grepl(tp53_loss2, BeatAML_aggrigate$Karyotype[i], fixed = T) & BeatAML_aggrigate$symbol[i] == "TP53") { BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]/2}
  
}

for(i in 1:nrow(BeatAML_aggrigate)){
  if(is.na(BeatAML_aggrigate$VAF_CN_corrected[i])) { BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$t_vaf[i]}
}

# correct for X-linked mutations in females
for(i in 1:nrow(BeatAML_aggrigate)){
  if(BeatAML_aggrigate$consensus_sex[i] == "Male" & BeatAML_aggrigate$chrom[i] == "X"){
    BeatAML_aggrigate$VAF_CN_corrected[i] = BeatAML_aggrigate$VAF_CN_corrected[i]/2
  }
}

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(BeatAML_aggrigate)){
  if(BeatAML_aggrigate$symbol[i] == "FLT3" & BeatAML_aggrigate$variant_class[i] == "ITD"){
    BeatAML_aggrigate$symbol[i] <- "FLT3-ITD"
  }
  if(BeatAML_aggrigate$symbol[i] == "FLT3" & BeatAML_aggrigate$variant_class[i] == "SNV"){
    BeatAML_aggrigate$symbol[i] <- "FLT3-TKD"
  }
  if(BeatAML_aggrigate$symbol[i] == "FLT3" & BeatAML_aggrigate$variant_class[i] == "deletion"){
    BeatAML_aggrigate$symbol[i] <- "FLT3-TKD"
  }
  if(BeatAML_aggrigate$symbol[i] == "FLT3" & BeatAML_aggrigate$variant_class[i] == "insertion"){
    BeatAML_aggrigate$symbol[i] <- "FLT3-ITD"
    BeatAML_aggrigate$variant_class[i] <- "ITD"
  }
  if(BeatAML_aggrigate$symbol[i] == "FLT3" & BeatAML_aggrigate$variant_class[i] == "tandem_duplication"){
    BeatAML_aggrigate$symbol[i] <- "FLT3-ITD"
    BeatAML_aggrigate$variant_class[i] <- "ITD"
  }
}


# VAF distribution ####

# subset to patients with required data types
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# de novo 
pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" & BeatAML_sample_data_types$isTransformed == "TRUE")

# pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE")

## filter to variants present in specific disease settings
BeatAML_aggrigate <- setDT(BeatAML_aggrigate)[LabId %chin% pt_subset$LabId]

mut_table <- aggregate(data.frame(count = BeatAML_aggrigate), list(value = BeatAML_aggrigate$symbol), length)
mut_table <- select(mut_table, "value", "count.symbol")
mut_table <- subset(mut_table, mut_table$count.symbol >= 5)

Beat_AML_muts <- setDT(BeatAML_aggrigate)[symbol %chin% mut_table$value] 

Beat_AML_muts[] <- lapply(Beat_AML_muts, gsub, pattern = "deletion", replacement = "Deletion", fixed = TRUE)
Beat_AML_muts[] <- lapply(Beat_AML_muts, gsub, pattern = "insertion", replacement = "Insertion", fixed = TRUE)
Beat_AML_muts[] <- lapply(Beat_AML_muts, gsub, pattern = "tandem_duplication", replacement = "ITD", fixed = TRUE)

Beat_AML_muts$VAF_CN_corrected=as.numeric(Beat_AML_muts$VAF_CN_corrected)

Beat_AML_muts$symbol <- with(Beat_AML_muts, reorder(symbol, -VAF_CN_corrected, median))


p = ggplot(Beat_AML_muts, aes(x=symbol, y=VAF_CN_corrected)) + 
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
  geom_jitter(aes(fill = variant_class), color = "black", shape = 21, position=position_jitter(0.2), size = 2.5) +
  # scale_fill_manual(values = c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")) +
  scale_fill_manual(values = c(Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF")) +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label = "VAF") +
  xlab(label = NULL) +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))
ggpar(p, legend.title = "Variant type")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/vaf_distribution_de_novo.pdf", dpi = 300, width = 8, height = 5, units = "in")


counts = as.data.frame(table(threshold_list_all$Gene))
names(counts)[1] = "Gene"

p =  ggplot(Beat_AML_muts, aes(x = symbol, y = VAF_CN_corrected)) +
  # ggtitle("Main Plot Title") +
  ylab("VAF") +
  xlab(NULL) +
  theme_cowplot(font_size = 15) +
  scale_shape_identity() +
  scale_fill_manual(values = c(Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF")) +
  geom_jitter(aes(fill = variant_class), color = "black", shape = 21, position=position_jitter(0.1), size = 2) +
  geom_flat_violin(position = position_nudge(x = 0.3, y = 0),
                   color = "grey", fill = "lightgrey",
                   adjust = 2,
                   alpha = 0.6, 
                   trim = TRUE, 
                   scale = "width") +
  geom_boxplot(position = position_nudge(x = 0.3, y = 0),
               notch = FALSE, 
               width = 0.2, 
               varwidth = FALSE, 
               outlier.shape = NA, 
               alpha = 0.3, 
               colour = "black", 
               show.legend = FALSE) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))  
ggpar(p, legend.title = "Variant type")

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/vaf_distribution_de_novo_raincloud.pdf", dpi = 300, width = 8, height = 5.5, units = "in")
write_csv(Beat_AML_muts, "~/Desktop/MetaAML_results/Figure_6/vaf_distribution_de_novo_raincloud.csv")





# plot the difference in VAF distribution per gene after CNA correction
raw_vaf = select(BeatAML_aggrigate, symbol, t_vaf)
raw_vaf$group = "VAF"


CNA_vaf = select(BeatAML_aggrigate, symbol, VAF_CN_corrected)
names(CNA_vaf)[2] = "t_vaf"
CNA_vaf$group = "VAF + \nCytogenetics"

all = rbind(raw_vaf, CNA_vaf)

# select mutations 
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1","NRAS", "KRAS", "FLT3-TKD", "FLT3-ITD", "PTPN11")

all = subset(all, all$symbol %in% mut_list)

all$group <- factor(all$group , levels=c("VAF", "VAF + \nCytogenetics"))


# Box plot facetted by "gene"
p <- ggpaired(all, x = "group", y = "t_vaf",
              color = "group", palette = "jama", 
              line.color = "lightgray", line.size = 0.4,
              facet.by = "symbol", short.panel.labs = T, ncol = 13)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE, label.x = 1.5, label.y = 75) + xlab(NULL) + ylab("VAF") + theme(legend.title=element_blank())
ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/Supplimental/BeatAML_vaf_correction_comparison_facet.pdf", dpi = 300, width = 15, height = 5, units = "in")


all$group <- factor(all$group , levels=c("VAF", "VAF + \nCytogenetics"))

p <- ggboxplot(all, x = "symbol", y = "t_vaf",
               color = "group", 
               palette = "jama",
               add = "jitter",
               shape = 20)
p + stat_compare_means(aes(group = group), label = "p.signif") + xlab(NULL) +
  # ylim(0,105) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + theme(legend.title=element_blank())

ggsave(filename = "~/Desktop/MetaAML_results/Figure_6/Supplimental/BeatAML_vaf_correction_comparison.pdf", dpi = 300, width = 15, height = 5, units = "in")



## Analyse all mutation relationships to a drug by performing linear regression between  VAF vs AUC 
# data ####
# read in drug sensitivity data
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# read in drug response data
drug_data <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 10)

# they have one drug that is not formated like all the others
drug_data[] <- lapply(drug_data, gsub, pattern = "17-AAG (Tanespimycin)", replacement = "Tanespimycin", fixed = TRUE)


## pair patient samples with drug response data

# first, remove rows that are second-hit mutations in the same gene at a different VAF
# arrange by VAF in descending order so that we keep only the highest VAF from overlapping hits
BeatAML_aggrigate_sub <- BeatAML_aggrigate %>%
  group_by(LabId) %>%
  arrange(desc(VAF_CN_corrected))

# create new data frame with the columns needed to determine mutation overlaps 
overlaps <- BeatAML_aggrigate_sub %>%
  select(LabId, symbol)

# create a column in BeatAML where TRUE if the value is duplicated earlier in the data table (does not put TRUE for first instance of that duplicate)
BeatAML_aggrigate_sub$duplicates <- NA
BeatAML_aggrigate_sub$duplicates <- duplicated(overlaps)

# keep only those that are true in keep column
BeatAML_aggrigate_sub <- BeatAML_aggrigate_sub[which(!BeatAML_aggrigate_sub$duplicates),]


# filter to only patients with both mutation and drug screening data
colnames(BeatAML_aggrigate_sub)[1] <- "lab_id"
drug_data <- setDT(drug_data)[lab_id %chin% BeatAML_aggrigate_sub$lab_id] 

drug_mut_all <- dplyr::right_join(drug_data, BeatAML_aggrigate_sub, by = "lab_id")

# loop through drugs and create a matrix of all drug-mutation-vaf combinations
drug_list <- list()
for(i in 1:nrow(drug_mut_all)){
  # print(i)
  drugmut <- as.character(drug_mut_all[i,1])
  drugmut <- gsub("\\s.*","",drugmut)
  drug_list[[i]] <- drugmut
}
drug_mut_list <- as.data.frame(do.call(rbind, drug_list))
names(drug_mut_list) <- "Inhibitor"

drug_mut_all <- cbind(drug_mut_list, drug_mut_all)
drug_mut_all$inhibitor <- NULL

drug_mut_all$auc <- as.numeric(drug_mut_all$auc)
drug_mut_all$VAF_CN_corrected <- as.numeric(drug_mut_all$VAF_CN_corrected)

mutations <- as.data.frame(unique(drug_mut_all$symbol))
inhibitors_list <- na.omit(as.data.frame(unique(drug_mut_all$Inhibitor)))


#### VAF - AUC regression ####
## Find cases where the vaf and auc correlate 
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/all")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/all/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/all/Sensitive")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/secondary")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/secondary/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/secondary/Sensitive")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/de_novo")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/de_novo/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/de_novo/Sensitive")
# decide on sample subset to run the analysis on
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# write a function that allows drug response analysis between VAF and drug AUC based on the sample subset of interest
drug_vaf_regression_function = function(sample_subset){
  
  # subset to patients with required data types
  if(sample_subset == "de_novo") {pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" &  BeatAML_sample_data_types$isDenovo == "TRUE" )}
  
  if(sample_subset == "secondary") {pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" &  BeatAML_sample_data_types$isTransformed == "TRUE" )}
  
  if(sample_subset == "all") {pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE")}
  
  ## filter to variants present in specific disease settings
  drug_mut <- setDT(drug_mut_all)[lab_id %chin% pt_subset$LabId]
  
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
    
    drug_sig <- as.character(inhibitors_list[i,1])
    for(j in 1:nrow(mutations)){
      # print(j)
      mut_sig <- as.character(mutations[j,1])
      
      drug_mut_sub <- as.data.frame(drug_mut[which(drug_mut$Inhibitor == drug_sig),])
      drug_mut_sub <- as.data.frame(drug_mut_sub[which(drug_mut_sub$symbol == mut_sig),])
      
      num_obs <- as.numeric(nrow(drug_mut_sub))
      
      if(num_obs >= 5){
        
        # find the delta in AUN and VAF for each case in order to filter later
        delta_vaf <- diff(range(drug_mut_sub$VAF_CN_corrected))
        min_vaf <- min(drug_mut_sub$VAF_CN_corrected)
        max_vaf <- max(drug_mut_sub$VAF_CN_corrected)
        delta_auc <- diff(range(drug_mut_sub$auc))
        min_auc <- min(drug_mut_sub$auc)
        max_auc <- max(drug_mut_sub$auc)
        
        # store all of the relationships in a dataframe
        if(num_obs >= 5 & delta_vaf >= 0.20 & delta_auc >= 50){
          
          # fit the data to a linear regression model
          lmf <- lm(drug_mut_sub$VAF_CN_corrected ~ drug_mut_sub$auc, data=drug_mut_sub)
          
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
          if(lmr >= 0.5 && lmp <= 0.05 && slp > 0 && delta_vaf >= 20 && delta_auc >= 50){
            
            # scatterplots ####
            p1 = ggscatter(drug_mut_sub, x = "VAF_CN_corrected", y = "auc",
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
                label.x = 10,
                label.y = (max_auc + 50)
              )
            
            vaf_auc_resistant_plots[[k]] = p1
            
            k = k + 1
            
            ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/",sample_subset,"/Resistant/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
          }
          if(lmr >= 0.5 && lmp <= 0.05 && slp < 0 && delta_vaf >= 20 && delta_auc >= 50){
            
            # make the scatterplot
            p2 = ggscatter(drug_mut_sub, x = "VAF_CN_corrected", y = "auc",
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
                label.x = 10,
                label.y = (max_auc + 50)
              )
            
            vaf_auc_sensitive_plots[[l]] = p2
            l = l + 1
            
            ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/",sample_subset,  "/Sensitive/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
          }
        }
      }
    }
  }
  beatAML_auc_vaf_comparison <- do.call(rbind, auc_vaf_list)
  
  # correct for multiple hypotheses by variable
  variables = unique(beatAML_auc_vaf_comparison$Inhibitor)
  
  var2_adj = list()
  a = 1
  
  for(i in 1:length(variables)){
    print(i)
    variable = as.character(variables[i])
    
    var2= subset(beatAML_auc_vaf_comparison, beatAML_auc_vaf_comparison$Inhibitor == variable)
    
    var2$q_val = p.adjust(var2$p_value)
    
    var2_adj[[a]] = var2 
    
    a = a + 1
    
  }
  beatAML_auc_vaf_comparison <- do.call(rbind, var2_adj)
  
  # beatAML_auc_vaf_comparison$fdr_corrected <- p.adjust(beatAML_auc_vaf_comparison$p_value, method = "fdr")
  # beatAML_auc_vaf_comparison$z_score <- qnorm(beatAML_auc_vaf_comparison$p_value)
  
  # write out results file
  write.csv(beatAML_auc_vaf_comparison, paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/", sample_subset,"/", sample_subset, "_drug_vaf_correlation.csv", sep = ""), row.names=FALSE)
  
  
  # make supplimental summary plots of all significant regressions
  # sensitive
  sensitive = list.clean(vaf_auc_sensitive_plots)
  
  p = cowplot::plot_grid(plotlist = sensitive, 
                         ncol = 7, align = "hv",   
                         rel_heights = c(1,1),
                         rel_widths = c(1,1))
  
  ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/", sample_subset,"/", sample_subset, "_sensitive_matrix.png", sep = ""), dpi = 150, width = 20, height = 5, units = "in")
  
  # resistant
  resistant = list.clean(vaf_auc_resistant_plots)
  
  p = cowplot::plot_grid(plotlist = resistant, 
                         ncol = 7, align = "hv",   
                         rel_heights = c(1,1),
                         rel_widths = c(1,1))
  
  ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/",sample_subset, "/", sample_subset, "_resistant_matrix.png", sep = ""), dpi = 150, width = 20, height = 7, units = "in")
  
  # summary plots ####
  # plot heatmaps of the delta AUC and delta VAF for the significant cases
  vaf_data <- subset(beatAML_auc_vaf_comparison, beatAML_auc_vaf_comparison$R_squared >= 0.25 & beatAML_auc_vaf_comparison$VAF_range >= 25 & beatAML_auc_vaf_comparison$AUC_range >= 75)
  
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
    if(vaf_data$q_val[i] < 0.1){
      vaf_data$star[i] = "*"
    }
  }
  
  # plot the heatmap for the most significant cases and manually correct a drug name
  vaf_data_sub = subset(vaf_data, vaf_data$p_value < 0.05)
  for(i in 1:nrow(vaf_data_sub)){
    if(vaf_data_sub$Inhibitor[i] == "Bay"){
      vaf_data_sub$Inhibitor[i] = "Bay 11-7085"
    }
  }
  
  p = ggplot(vaf_data_sub, aes(reorder(Mutated_Gene, -count.Inhibitor), reorder(Inhibitor, count.Inhibitor), colour = AUC_range, size = VAF_range, label = star)) +
    geom_point() +
    theme_cowplot() +
    scale_color_gradient2(low = "#b2182b", high = "#2166ac", mid = "#f7f7f7",
                          name=expression(Delta~"AUC")) +
    geom_point(shape = 21, color = "black") +
    # scale_size_area(max_size = 5,breaks=c(25,50,75)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "#f0f0f0"))+
    xlab(label= NULL) +
    ylab(label="Inhibitor") +
    labs(title = NULL) +
    geom_text(size = 7.5, color = "#525252",  hjust = 0, nudge_x = 0.25)+
    theme(legend.position=c(0.8,0.85)) 
  
  g = guide_legend(override.aes=list(colour="lightgrey"), expression(Delta~"VAF"))
  p + guides(size = g) +
    annotate(geom="text", x = 10.5, y = 27, label="* q < 0.1",
             color="#525252")
  
  ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/", sample_subset, "/", sample_subset, "_drug_vaf_correlation_summary_plot.pdf", sep = ""), dpi = 300, width = 5, height =  11, units = "in")
  
  
  # plot the heatmap for everything
  p = ggplot(vaf_data, aes(reorder(Mutated_Gene, -count.Inhibitor), reorder(Inhibitor, count.Inhibitor), colour = AUC_range, size = VAF_range, label = star)) +
    geom_point() +
    theme_cowplot() +
    scale_color_gradient2(low = "#b2182b", high = "#2166ac", mid = "#f7f7f7",
                          name=expression(Delta~"AUC")) +
    geom_point(shape = 21, color = "black") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "#f0f0f0"))+
    xlab(label= NULL) +
    ylab(label="Inhibitor") +
    labs(title = NULL) +
    geom_text(size = 7.5, color = "#525252",  hjust = 0, nudge_x = 0.25) +
    theme(legend.position=c(0.95,0.65)) 
  
  g = guide_legend(override.aes=list(colour="lightgrey"), expression(Delta~"VAF"))
  p + guides(size = g) +
    annotate(geom="text", x = 3, y = 85, label="* q < 0.1",
             color="#525252") +
    coord_flip()
  
  ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/", sample_subset, "/", sample_subset, "_drug_vaf_correlation_all_summary_plot.pdf", sep = ""), dpi = 300, width = 20, height =  5, units = "in")
}

drug_vaf_regression_function("de_novo")
drug_vaf_regression_function("secondary")
drug_vaf_regression_function("all")

# binary analysis ####
## loop through all gene/drug interactions and calculate the p-value for a t.test between mutated and non-mutated samples for each drug
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/all")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/all/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/all/Sensitive")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/secondary")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/secondary/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/secondary/Sensitive")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/de_novo")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/de_novo/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/de_novo/Sensitive")

# load patients with avaliable data type annotations
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# create a dataframe of all the unique drugs used in the screen
drug_list <- na.omit(as.data.frame(unique(drug_mut_all$Inhibitor)))
# 122 unique drugs

drug_mutation_binary_function = function(sample_subset, individual_plots){
  
  # subset to patients with required data types
  if(sample_subset == "de_novo") {pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" &  BeatAML_sample_data_types$isDenovo == "TRUE" & BeatAML_sample_data_types$totalDrug == "y")}
  
  if(sample_subset == "secondary") {pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" &  BeatAML_sample_data_types$isTransformed == "TRUE"  & BeatAML_sample_data_types$totalDrug == "y")}
  
  if(sample_subset == "all") {pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" & BeatAML_sample_data_types$totalDrug == "y")}
  
  ## filter to variants present in specific disease settings
  drug_mut <- setDT(drug_mut_all)[lab_id %chin% pt_subset$LabId]
  
  datalist_beatAML <- list()
  z <- 1
  for(i in 1:nrow(drug_list)){
    print(i)
    
    drug <- as.character(drug_list[i,1])
    subdat <- as.data.frame(subset(drug_mut, drug_mut$Inhibitor == drug))
    
    # filter to unique patients
    subdat_uniq <- subset(subdat, !duplicated(lab_id))
    
    # mut_table frome line 198
    for(j in 1:nrow(mut_table)){
      # print(j)
      gene_name <- as.character(mut_table[j,1])
      
      # create column to deliniate if the patient has a somatic mutation of interest
      subdat_uniq2 = subdat_uniq %>%
        select(Inhibitor, lab_id, auc)
      
      # find pts with a mutation in the gene of interest
      pt_mut = setDT(drug_mut)[symbol %chin% mut_table[j,1]]
      pt_mut = pt_mut %>% 
        select(lab_id, ) %>%
        unique()
      
      if(nrow(pt_mut) >= 5){
        subdat_uniq2$mutated = ifelse(subdat_uniq2$lab_id %chin% pt_mut$lab_id, gene_name, "WT")
        
        # calculate the mean drug AUC for the mutated and WT samples
        mut <- subdat_uniq2[which(subdat_uniq2$mutated == gene_name),]
        nummut <- nrow(mut)
        mut <-  as.numeric(mean(mut$auc))
        
        nmut <- subdat_uniq2[which(subdat_uniq2$mutated == "WT"),]
        nmut <-  as.numeric(mean(nmut$auc))
        
        delta <- as.numeric(mut - nmut)
        
        # if there are enough cases for analysis, calculate differences in AUC between mutated and WT groups, taking into account the distribution of the data to inform the statistical test used
        
        if(nummut >= 5){
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
        # plot boxplot for distribution betweeen mut and wt samples for the drug in the loop
        if(individual_plots == "yes"){
          
          if(p_value < 0.05){
            subdat_uniq2$mutated = ifelse(subdat_uniq2$mutated != "WT", "Mut", "WT")
            
            if(delta > 0){
              
              subdat_uniq2$mutated <- factor(subdat_uniq2$mutated , levels=c("WT", "Mut"))
              
              ggplot(subdat_uniq2, aes(x = mutated, y = auc)) +
                ylab(label= "AUC") +
                xlab(label = NULL) +
                theme_cowplot(font_size = 12) +
                stat_compare_means() +
                scale_fill_manual(values = c("WT"="lightgrey", "Mut" = "#b2182b")) +
                labs(title = paste(drug, " vs. ", gene_name, sep = "")) +
                ylim((min(subdat_uniq2$auc)), (max(subdat_uniq2$auc) + 15)) +
                theme(legend.position="none")  +
                geom_point(aes(fill = mutated), color = "black", shape = 21, position=position_jitter(0.15), size = 2) +
                geom_flat_violin(position = position_nudge(x = 0.3, y = 0),
                                 color = "grey", fill = "lightgrey",
                                 adjust = 2,
                                 alpha = 0.6, 
                                 trim = TRUE, 
                                 scale = "width") +
                geom_boxplot(position = position_nudge(x = 0.3, y = 0),
                             notch = FALSE, 
                             width = 0.2, 
                             varwidth = FALSE, 
                             outlier.shape = NA, 
                             alpha = 0.3, 
                             colour = "black", 
                             show.legend = FALSE)
              
              ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/",sample_subset,"/Sensitive/",drug,"_",gene_name,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
            }
            if(delta < 0){
              
              subdat_uniq2$mutated <- factor(subdat_uniq2$mutated , levels=c("WT", "Mut"))
              
              ggplot(subdat_uniq2, aes(x = mutated, y = auc)) +
                ylab(label= "AUC") +
                xlab(label = NULL) +
                theme_cowplot(font_size = 12) +
                stat_compare_means() +
                scale_fill_manual(values = c("WT"="lightgrey", "Mut" = "#4393c3")) +
                labs(title = paste(drug, " vs. ", gene_name, sep = "")) +
                ylim((min(subdat_uniq2$auc)), (max(subdat_uniq2$auc) + 15)) +
                theme(legend.position="none")  +
                geom_point(aes(fill = mutated), color = "black", shape = 21, position=position_jitter(0.15), size = 2) +
                geom_flat_violin(position = position_nudge(x = 0.3, y = 0),
                                 color = "grey", fill = "lightgrey",
                                 adjust = 2,
                                 alpha = 0.6, 
                                 trim = TRUE, 
                                 scale = "width") +
                geom_boxplot(position = position_nudge(x = 0.3, y = 0),
                             notch = FALSE, 
                             width = 0.2, 
                             varwidth = FALSE, 
                             outlier.shape = NA, 
                             alpha = 0.3, 
                             colour = "black", 
                             show.legend = FALSE)

              ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/",sample_subset,"/Resistant/",drug,"_",gene_name,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
            }
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
      beatAML_sensitivity_difference$color[i] <- "ns."
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
  
  write.csv(beatAML_sensitivity_difference, paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/",sample_subset,"/binary_drug_results.csv", sep = ""), row.names=FALSE)
  
  # beatAML_sensitivity_difference = read.csv("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/de_novo/binary_drug_results.csv")
  
  ## plot distribution results as a volcano plot
  p = ggplot(beatAML_sensitivity_difference, aes(x=auc_diff, y=-log10(fdr_adjusted), color = factor(color), size = count.symbol)) +
    theme_cowplot() +
    geom_point(alpha = 0.75) + 
    geom_point(shape = 21, color = "darkgrey", alpha = 0.25) +
    geom_label_repel(aes(label=point_label), size = 3,  show_guide = F) +
    scale_colour_manual(values = c("Sensitive"= "#b2182b", "Resistant"="#2166ac", "ns."="grey")) + 
    theme(legend.position="right") +
    geom_hline(yintercept = 1.30103,  linetype = "dashed", color = "grey") +
    ylab(label= "-log10(FDR) q-value") +
    xlab(label= expression(Delta~"AUC (wt - mut)")) +
    labs(title = NULL) +
    theme(plot.title = element_text(color="black", size=20)) 
  
  g = guide_legend(override.aes=list(colour="grey"), "Number of\nsamples")
  h = guide_legend("")
  
  p + guides(size = g, color = h) +
    annotate(geom="text", x = 100, y = 1.5, label="q < 0.1",
             color="grey")
  
  ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/",sample_subset,"/binary_drug_results.pdf", sep = ""), dpi = 300, width = 7, height = 4.5, units = "in")
  
}

drug_mutation_binary_function(sample_subset = "de_novo", individual_plots = "yes")
drug_mutation_binary_function("secondary", individual_plots = "yes")
drug_mutation_binary_function("all", individual_plots = "yes")


# delete large dataframes
rm(drug_list)
rm(datalist_beatAML)



# overlap ####
# overlap of significant hits for both approaches using an upset plot
binary_results=read.csv("~/Desktop/MetaAML_results/Figure_6/Supplimental/drug_vaf_correlation/binary/de_novo/binary_drug_results.csv")
binary_results=subset(binary_results, binary_results$fdr_adjusted <= 0.1)
binary_sensitive=subset(binary_results, binary_results$auc_diff < 0)
binary_resistant=subset(binary_results, binary_results$auc_diff > 0)


vaf_results=read.csv("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation.csv")
vaf_results=subset(vaf_results, vaf_results$p_value <= 0.05 & vaf_results$R_squared >= 0.5 & vaf_results$VAF_range >= 0.25 & vaf_results$AUC_range >= 75)
vaf_sensitive=subset(vaf_results, vaf_results$Slope < 0)
vaf_resistant=subset(vaf_results, vaf_results$Slope > 0)

# find overlap for interactions
b_v_s_overlap=inner_join(binary_sensitive, vaf_sensitive, by = c("Inhibitor", "Mutated_Gene"))
b_v_r_overlap=inner_join(binary_resistant, vaf_resistant, by = c("Inhibitor", "Mutated_Gene"))


upset <- data.frame(matrix(0, nrow = 100, ncol = 4))
names(upset) <- c("binary_sensitive", "binary_resistant", "vaf_sensitive", "vaf_resistant")

for(i in 1:nrow(upset)){
  if(i <= 37){
    upset$binary_sensitive[i] = 1
  }
  if(i > 37 & i <= 57){
    upset$binary_resistant[i] = 1
  }
  if(i > 57 & i <= 73){
    upset$vaf_sensitive[i] = 1
  }
  if(i > 73){
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



# co-occuring mutations and FLT3 vaf regression ####
dir.create("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/FLT3_specific")
dir.create("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/FLT3_specific/Resistant")
dir.create("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/FLT3_specific/Sensitive")


pt_subset <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$isRelapse == "FALSE" &  BeatAML_sample_data_types$isDenovo == "TRUE" & BeatAML_sample_data_types$totalDrug == "y")

drug_mut <- setDT(drug_mut_all)[lab_id %chin% pt_subset$LabId]

# select the different flt3-mutated samples
flt3_itd_samples = subset(drug_mut, drug_mut$symbol == "FLT3-ITD")
flt3_tkd_samples = subset(drug_mut, drug_mut$symbol == "FLT3-TKD")


background_genotype = c("DNMT3A", "NMP1", "NRAS", "TET2", "RUNX1")

auc_vaf_list <- list()
vaf_auc_sensitive_plots = list()

z = 1
k = 1

for(i in 1:length(background_genotype)){
  print(background_genotype[i])
  
  # subset to only include patients with co-occuring mutations in the background genotype
  drug_mut_sub = subset(drug_mut, drug_mut$symbol == background_genotype[i]) 
  
  # subset based on flt3-itd cases
  drug_mut_sub = setDT(drug_mut_sub)[lab_id %chin% flt3_itd_samples$lab_id]
  drug_mut_sub = select(drug_mut_sub, lab_id, symbol, VAF_CN_corrected)
  names(drug_mut_sub) = c("lab_id", "gene2", "gene2_vaf")
  
  itd = inner_join(flt3_itd_samples, drug_mut_sub, by = "lab_id")
  
  # subset based on flt3-tkd cases
  drug_mut_sub = subset(drug_mut, drug_mut$symbol == background_genotype[i]) 
  drug_mut_sub = setDT(drug_mut_sub)[lab_id %chin% flt3_tkd_samples$lab_id]
  drug_mut_sub = select(drug_mut_sub, lab_id, symbol, VAF_CN_corrected)
  names(drug_mut_sub) = c("lab_id", "gene2", "gene2_vaf")
  
  tkd = inner_join(flt3_tkd_samples, drug_mut_sub, by = "lab_id")
  
  data_list = list(itd,tkd)
  
  # auc_vaf_list <- list()
  # vaf_auc_sensitive_plots = list()
  # 
  # z = 1
  # k = 1
  
  for(a in 1:length(data_list)){
    
    sub_data = as.data.frame(data_list[a])
    flt3_type = sub_data$symbol[1]
    
    for(j in 1:nrow(inhibitors_list)){
      drug_sig <- as.character(inhibitors_list[j,1])
      
      sub = unique(subset(sub_data, sub_data$Inhibitor == drug_sig))
      
      sub$auc = as.numeric(sub$auc)
      
      mut_sig = sub$gene2[1]
      
      if(nrow(sub) >= 5){
        # print(i)
        # find the delta in AUN and VAF for each case in order to filter later
        delta_vaf <- as.numeric(diff(range(sub$VAF_CN_corrected)))
        min_vaf <- as.numeric(min(sub$VAF_CN_corrected))
        max_vaf <- as.numeric(max(sub$VAF_CN_corrected))
        delta_auc <- as.numeric(diff(range(sub$auc)))
        min_auc <- as.numeric(min(sub$auc))
        max_auc <- as.numeric(max(sub$auc))
        
        # store all of the relationships in a dataframe
        if(delta_vaf >= 0.2 & delta_auc >= 50){
          
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
          if(lmr >= 0.5 && lmp <= 0.05 && slp > 0 && delta_vaf >= 0.2 && delta_auc >= 50){
            
            # scatterplots ####
            p1 = ggscatter(sub, x = "VAF_CN_corrected", y = "auc",
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
            
            ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/FLT3_specific/Resistant/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
          }
          if(lmr >= 0.5 && lmp <= 0.05 && slp < 0 && delta_vaf >= 0.2 && delta_auc >= 50){
            #   
            # make the scatterplot
            p2 = ggscatter(sub, x = "VAF_CN_corrected", y = "auc",
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
            
            ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/FLT3_specific/Sensitive/",drug_sig,"_",mut_sig,".pdf", sep = ""), dpi = 300, width = 3, height = 3, units = "in")
          }
        }
      }
    }  
  }
  
}
# 
flt3_auc_vaf_comparison <- do.call(rbind, auc_vaf_list)
flt3_auc_vaf_comparison$fdr_corrected <- p.adjust(flt3_auc_vaf_comparison$p_value, method = "fdr")


# write out results file
write.csv(flt3_auc_vaf_comparison, "~/Desktop/MetaAML_results/Figure_6/drug_vaf_correlation/FLT3_specific/flt3_drug_vaf_correlation.csv", row.names=FALSE)
