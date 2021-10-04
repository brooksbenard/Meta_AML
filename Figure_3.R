# ========================================================================================================================================== #
# Figure_3.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 08/23/2021
# Description: This script will perform survival analyses based on VAF thresholds as seen in Figure 3 of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# ================ #
# Load packages ####
# ================ #
# Package names
packages <- c("ggplot2",  "dplyr", "cometExactTest", "stringr", "maditr", "reshape2", "data.table", "epitools", "corrplot", "plyr", "muhaz", "reshape", "survival", "survivalAnalysis", "survMisc", "survminer", "ggsci", "vegan", "ggrepel", "ggforce", "rstatix", "effsize", "psych", "maxstat", "RCurl", "ggpubr", "UpSetR", "cowplot", "readxl", "scales", "rlist", "tidyverse", "PupillometryR")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# create directories
dir.create("~/Desktop/MetaAML_results/Figure_3")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental")

#### CNA correction ####
# Using the cytogenetic data, recalculate the VAF for mutations falling on focal or broad CNAs
# Pappamanuel/Gerstung VAF correction is performed using code provided by the authors
# TCGA VAFs are corrected using a gene-level CNA supplimentary file
# for all other cohorts, manual correction is performed using reported cytogenetic/karyotypic data. 

#### Papaemmanuil/Gerstung ####
# download.file("https://github.com/gerstung-lab/AML-multistage/blob/master/data/AMLSG_Clinical_Anon.RData?raw=true", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")

load("~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")

# 1.3.2.2 Mutation data
# download.file("https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/AMLSG_Genetic.txt", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Genetic.txt")

mutationData = read.table("~/Desktop/MetaAML_results/raw_data/AMLSG_Genetic.txt", sep="\t", header=TRUE, strip.white = TRUE) 
mutationData$SAMPLE_NAME <- factor(as.character(mutationData$SAMPLE_NAME), levels = levels(clinicalData$PDID)) 
## Refactor
mutationTable <- (table(mutationData[mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & mutationData$FINAL_CALL == "OK" ,c("SAMPLE_NAME","GENE")]) > 0)+0
dim(mutationTable)

all(rownames(mutationTable)==clinicalData$PDID)

# 1.3.4 Covariates
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
                Cytogenetics = clinicalData[,grep("^(t_)|(inv)|(abn)|(plus)|(minus)|(mono)|(complex)",colnames(clinicalData))],
                Nuisance = data.frame( as.integer(clinicalData$Study), 
                                       Date=scale(as.numeric(clinicalData$ERDate), scale=FALSE), 
                                       MissingCyto=is.na(clinicalData$t_15_17)+0),
                Demographics = clinicalData[,c("AOD","gender")],
                Clinical = cbind(clinicalData[, c("Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH","HB","platelet",
                                                  "Splenomegaly")], as.integer(clinicalData$TypeAML)),
                MolRisk = as.integer(clinicalData$M_Risk))
dataList$Genetics$CEBPA <- clinicalData$CEBPA # encoded as 0,1,2 
dataList$Genetics$CEBPA_mono <- clinicalData$CEBPA == 1 # encoded as 0,1,2 
dataList$Genetics$CEBPA_bi <- clinicalData$CEBPA == 2 # encoded as 0,1,2 
dataList$Genetics$CEBPA <- NULL
dataList$Genetics$FLT3 <- NULL
dataList$Genetics$FLT3_ITD <- clinicalData$FLT3_ITD != "0"
dataList$Genetics$FLT3_TKD <- clinicalData$FLT3_TKD != "0"
dataList$Genetics$FLT3_other <- clinicalData$FLT3_other != "0"
dataList$Genetics$IDH2_p172 <- table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("172", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2_p140 <- table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("140", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2 <- NULL
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Cytogenetics$MLL_PTD <- NULL
dataList$Genetics = dataList$Genetics + 0

# Condensing to a data.frame
dataRaw <- do.call(cbind,dataList)
names(dataRaw) <- unlist(sapply(dataList, names)) 
dataFrame = dataRaw
dim(dataFrame)
rownames(dataFrame) <- clinicalData$PDID

# download.file("https://github.com/gerstung-lab/AML-multistage/raw/master/data/AMLSG_Karyotypes.txt", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Karyotypes.txt")

# 1.3.5 Subclonal mutations
copyNumbers = cbind(dataList$Cytogenetics[grep(c("minus|plus|mono"), colnames(dataList$Cytogenetics))], clinicalData$gender)
copyNumbers$minus7 <- (copyNumbers$minus7 | copyNumbers$minus7q) + 0
copyNumbers$minus7q <- NULL

for(i in 1:ncol(copyNumbers)){
  if(grepl("plus", colnames(copyNumbers)[i])){
    copyNumbers[,i] = copyNumbers[,i]*3
  }
}
copyNumbers[copyNumbers==0 | is.na(copyNumbers)] = 2
colnames(copyNumbers) = c(5,7,8,9,12,13,17,18,20,21,22,"Y",11,4,"X")
rownames(copyNumbers) <- clinicalData$PDID

minusY = dataList$Cytogenetics$minusY[is.na(dataList$Cytogenetics$minusY)] <- 0

copyNumbers$Y <- c(1:0)[clinicalData$gender] - minusY
cn = sapply(1:nrow(mutationData), function(i) {
  c=copyNumbers[mutationData$SAMPLE_NAME[i],match(mutationData$CHR[i ], colnames(copyNumbers))]; if(length(c)==0) 2 else c
}
)
vaf <- as.numeric(as.character(mutationData$X._MUT_IN_TUM))

depth <- as.numeric(as.character(mutationData$TUM_DEPTH))


CNA = data.frame(cn)
rownames(CNA) = rownames(mutationData)
colnames(CNA) = "CNA"

mutationData$CN_status = as.numeric(CNA$CNA)
mutationData$VAF = as.numeric(mutationData$VAF)
# mutationData$CN_status =mutationData$CN_status %>% replace_na(2)
mutationData$CN_status[is.na(mutationData$CN_status)] <- 2
mutationData$VAF_CN_corrected = mutationData$VAF

for(i in 1:nrow(mutationData)){
  if(!is.na(mutationData$VAF_CN_corrected[i])){
  if(mutationData$CN_status[i] == 1){
    mutationData$VAF_CN_corrected[i] = mutationData$VAF[i]/2
  }
  if(mutationData$CN_status[i] == 2 & mutationData$VAF[i] > 55){
    mutationData$VAF_CN_corrected[i] = mutationData$VAF[i]/2
  }
  if(mutationData$CN_status[i] == 3 & mutationData$VAF[i] <= 35){
    mutationData$VAF_CN_corrected[i] = mutationData$VAF[i]*1.5
  }
  if(mutationData$CN_status[i] == 3 & mutationData$VAF[i] > 35){
    mutationData$VAF_CN_corrected[i] = mutationData$VAF[i]*0.75
  }
  }
}

for(i in 1:nrow(mutationData)){
  if(!is.na(mutationData$VAF_CN_corrected[i])){
    if(mutationData$VAF_CN_corrected[i] > 100){
      mutationData$VAF_CN_corrected[i] = 100
    }
  }
}

papaemmanuil_muts = mutationData

# for some reason, this cohort uses SFRS2 instead of SRSF2 like all the other cohorts
for(i in 1:nrow(papaemmanuil_muts)){
  if(papaemmanuil_muts$GENE[i] == "SFRS2"){
    papaemmanuil_muts$GENE[i] <- "SRSF2"
  }
}

papaemmanuil_muts = select(papaemmanuil_muts, SAMPLE_NAME, GENE, VAF, VAF_CN_corrected)
names(papaemmanuil_muts) = c("Sample", "Gene", "VAF", "VAF_CN_corrected")

papaemmanuil_muts$VAF <- round(papaemmanuil_muts$VAF, 2)

load("~/Desktop/MetaAML_results/final_data_matrix.RData")
papaemmanuil_sub = subset(final_data_matrix, final_data_matrix$Cohort == "Papaemmanuil")
papaemmanuil_corrected = left_join(papaemmanuil_sub, papaemmanuil_muts, by=c("Sample","Gene", "VAF"))


#### TCGA ####
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
tcga_dat = subset(final_data_matrix, final_data_matrix$Cohort == "TCGA")

# read in copy number alteration data frame
CNA <- as.data.frame(read.table("~/Desktop/MetaAML_results/raw_data/data_CNA.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE))
CNA$Entrez_Gene_Id = NULL

# reformat colnames to match those in the mutation dataframe
colnames(CNA) = gsub('\\.', '-', colnames(CNA))
colnames(CNA) = as.character(substr(colnames(CNA), 1, nchar(colnames(CNA))-3))

tcga_dat$cn = NA

for (i in 1:nrow(tcga_dat)) {
  current_id <- tcga_dat$Sample[i]
  current_symbol <- tcga_dat$Gene[i]
  delta <- as.integer(CNA[which(CNA$Hugo_Sym == current_symbol), which(colnames(CNA) == current_id)])
  if (length(delta) == 0) {
    delta <- 0
  }
  tcga_dat$cn[i] <- 2 + delta
}

tcga_dat$VAF_CN_corrected = NA


for(i in 1:nrow(tcga_dat)){
  if(!is.na(tcga_dat$VAF[i])){
    if(tcga_dat$cn[i] == 1){
      tcga_dat$VAF_CN_corrected[i] = tcga_dat$VAF[i]/2
    }
    if(tcga_dat$cn[i] == 2 & tcga_dat$VAF[i] <= 55){
      tcga_dat$VAF_CN_corrected[i] = tcga_dat$VAF[i]
    }
    if(tcga_dat$cn[i] == 2 & tcga_dat$VAF[i] > 55){
      tcga_dat$VAF_CN_corrected[i] = tcga_dat$VAF[i]/2
    }
    if(tcga_dat$cn[i] == 3 & tcga_dat$VAF[i] <= 55){
      tcga_dat$VAF_CN_corrected[i] = tcga_dat$VAF[i]*1.5
    }
    if(tcga_dat$cn[i] == 3 & tcga_dat$VAF[i] > 55){
      tcga_dat$VAF_CN_corrected[i] = tcga_dat$VAF[i]*0.75
    }
  }
}

for(i in 1:nrow(tcga_dat)){
  if(is.na(tcga_dat$VAF_CN_corrected[i])) { tcga_dat$VAF_CN_corrected[i] = tcga_dat$VAF[i]}
}


#### Other cohorts ####

# download a file correlating the chromosome loci with gene ID
# http://uswest.ensembl.org/biomart/martview/a3043f537a26692111e9a7a94003ff68
all_genes <- read.table("~/Downloads/mart_export.txt", sep = "\t", header = T, fill = T,  stringsAsFactors = FALSE, quote = "")

mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1","NRAS", "KRAS", "FLT3", "PTPN11")

all_genes = subset(all_genes, all_genes$Gene.name %in% mut_list)

all_genes$Karyotype_arm = gsub("\\..*","",all_genes$Karyotype.band)
all_genes$partial_annotation = paste("(",all_genes$Chromosome.scaffold.name,")","(",all_genes$Karyotype_arm,")", sep = "")
all_genes$full_annotation = paste("(",all_genes$Chromosome.scaffold.name,")","(",all_genes$Karyotype.band,")", sep = "")

# append to mutation file
names(all_genes)[1] = "Gene"

# load and select data from the aggregated cohort
# VAf distribution ####
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

vaf_sub = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
# vaf_sub = subset(vaf_sub, vaf_sub$Cohort == "Tyner")

vaf_sub = subset(vaf_sub, vaf_sub$Cohort != "TCGA" & vaf_sub$Cohort != "Papaemmanuil")
vaf_sub = subset(vaf_sub, vaf_sub$mut_freq_gene >= 50)
vaf_sub = subset(vaf_sub, vaf_sub$Gene != "MLL")

# append the chromosomal arm/band column to the mutation table
cohort_aggrigate = left_join(vaf_sub, all_genes, by = "Gene")

cohort_aggrigate$Cytogenetics = tolower(cohort_aggrigate$Cytogenetics)

# subset mutations to those where karyotyping data suggests deletions or amplifications
cohort_aggrigate$Cytogenetics = gsub('\\s+', '', cohort_aggrigate$Cytogenetics)
cohort_aggrigate$VAF_CN_corrected = NA

for(i in 1:nrow(cohort_aggrigate)){
  if(!is.na(cohort_aggrigate$VAF_male_x[i])){
    # find whole chromosome gains or losses
    ch_gain = paste("+",cohort_aggrigate$Chromosome.scaffold.name[i], ",", sep = "")
    ch_loss = paste("-",cohort_aggrigate$Chromosome.scaffold.name[i], ",", sep = "")
    
    # correct for VAF based on copy number gains for that chromosome  
    if(grepl(ch_gain, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$VAF_male_x[i] > 35) {cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF_male_x[i]*0.75}
    if(grepl(ch_gain, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$VAF_male_x[i] <= 35) {cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF_male_x[i]*1.5}
    # correct for VAF based on copy number loss for that chromosome
    if(grepl(ch_loss, cohort_aggrigate$Cytogenetics[i])){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF_male_x[i]/2}
    
    # find focal chromosome gains or losses
    locus_gain = paste("add",cohort_aggrigate$partial_annotation[i],"|","add",cohort_aggrigate$full_annotation[i], sep = "")
    locus_loss = paste("del",cohort_aggrigate$partial_annotation[i],"|","del",cohort_aggrigate$full_annotation[i], sep = "")
    
    # correct for VAF based on broad copy number gains at the gene locus  
    if(grepl(locus_gain, cohort_aggrigate$Cytogenetics[i], fixed = T) & cohort_aggrigate$VAF_male_x[i] > 35){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF_male_x[i]*0.75}
    if(grepl(locus_gain, cohort_aggrigate$Cytogenetics[i], fixed = T) & cohort_aggrigate$VAF_male_x[i] <= 35){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*1.5}
    # correct for VAF based on broad copy number loss at the gene locus
    if(grepl(locus_loss, cohort_aggrigate$Cytogenetics[i], fixed = T)){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAVAF_male_xF[i]/2}
  }
}

# TP53 is a unique case. Code manually for focal deletions
tp53_loss1 = "del(17)(p13)"
tp53_loss2 = "del(17)(p11.2p13)"

for(i in 1:nrow(cohort_aggrigate)){
  if(grepl(tp53_loss1, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$Gene[i] == "TP53") { cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2}
  if(grepl(tp53_loss2, cohort_aggrigate$Cytogenetics[i], fixed = T) & cohort_aggrigate$Gene[i] == "TP53") { cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2}
}
# if no copy number differenes detected, populate the raw VAF
for(i in 1:nrow(cohort_aggrigate)){
  if(is.na(cohort_aggrigate$VAF_CN_corrected[i])) { cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF_male_x[i]}
}

# bind all of the cohorts together into a final data frame
# remove the copy number column from some dataframes to make it compatible with the other data frames
tcga_dat$cn = NULL
cohort_aggrigate[,28:32] = NULL

# bind all cohorts together
final_data_matrix = rbind(papaemmanuil_corrected, tcga_dat, cohort_aggrigate, fill = TRUE)
final_data_matrix_2 = subset(final_data_matrix, final_data_matrix$Subset == "de_novo" & final_data_matrix$VAF_CN_corrected > 0)
save(final_data_matrix_2,  file = "~/Desktop/MetaAML_results/final_data_matrix_2.RData")

# delineate FLT3 mutations
for(i in 1:nrow(final_data_matrix_2)){
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] %in% c("ITD", "INDEL")){
    final_data_matrix_2$Gene[i] <- "FLT3-ITD"
  }
  if(final_data_matrix_2$Gene[i] == "FLT3" & final_data_matrix_2$variant_type[i] %in% c("SNV", "Deletion", "other")){
    final_data_matrix_2$Gene[i] <- "FLT3-TKD"
  }
}


# plot the difference in VAF distribution per gene after CNA correction
raw_vaf = select(final_data_matrix_2, Gene, VAF, Cohort)
raw_vaf$group = "VAF"

CNA_vaf = select(final_data_matrix_2, Gene, VAF_CN_corrected, Cohort)
names(CNA_vaf)[2] = "VAF"
CNA_vaf$group = "VAF + \nCytogenetics"

all = rbind(raw_vaf, CNA_vaf)

# select mutations 
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1", "FLT3-ITD","NRAS", "KRAS", "FLT3-TKD", "PTPN11")

all = subset(all, all$Gene %in% mut_list)

all$group = factor(all$group , levels=c("VAF", "VAF + \nCytogenetics"))

genes = sort(unique(all$Gene))

# Box plot facetted by "gene"
p <- ggpaired(all, x = "group", y = "VAF",
              color = "group", palette = "jama", 
              line.color = "lightgray", line.size = 0.4,
              facet.by = "Gene", short.panel.labs = T)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE) + xlab(NULL) + ylab("VAF") + theme(legend.title=element_blank())
ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/all_vaf_correction_comparison_facet_gene.pdf", dpi = 300, width = 10, height = 12, units = "in")

# plot total vaf distribution for recurrent mutations
p <- ggboxplot(all, x = "Gene", y = "VAF",
               color = "group", palette = "jama",
               add = "jitter",
               order = genes, 
               shape = 21)
p + stat_compare_means(aes(group = group), label = "p.signif", hide.ns = TRUE) + xlab(NULL) +
  ylim(0,110) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + theme(legend.title=element_blank())

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/all_vaf_correction_comparison.pdf", dpi = 300, width = 15, height = 3, units = "in")


# now plot the VAF distribution per gene colored by if the mutation is above or below the median VAF for that gene
genes = unique(final_data_matrix_2$Gene)

final_data_matrix_2$median_threshold = NA

threshold_list = list()
z = 1

for(i in 1:length(genes)){
  # select the gene of interest
  sub = subset(final_data_matrix_2, Gene == genes[i])
  # calculate the median vaf based on the X-corrected VAF
  med_vaf = as.numeric(median(sub$VAF_CN_corrected, na.rm = T))
  
  sub$median_threshold = ifelse(sub$VAF >= med_vaf, "High", "Low")
  sub$median_threshold_NA = ifelse(sub$VAF_CN_corrected >= med_vaf, "High", "Low")
  
  threshold_list[[i]] = data.frame(sub)
}

threshold_list_all = do.call(rbind, threshold_list)
rm(threshold_list)

threshold_list_all = threshold_list_all[!is.na(threshold_list_all$VAF_CN_corrected),]

# select mutations 
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1", "FLT3-ITD","NRAS", "KRAS", "FLT3-TKD", "PTPN11")

threshold_list_all = subset(threshold_list_all, threshold_list_all$Gene %in% mut_list)

threshold_list_all$Gene <- with(threshold_list_all, reorder(Gene, -VAF_CN_corrected, median))

counts = as.data.frame(table(threshold_list_all$Gene))
names(counts)[1] = "Gene"

p =  ggplot(threshold_list_all, aes(x = Gene, y = VAF_CN_corrected)) +
    # ggtitle("Main Plot Title") +
    ylab("VAF") +
    xlab(NULL) +
    theme_cowplot(font_size = 10) +
    scale_shape_identity() +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_manual(values = c("#cb181d", "#3690c0")) +
    geom_point(aes(fill = median_threshold_NA), alpha = 0.5, color = "black", shape = 21, position=position_jitter(0.15), size = 1) +
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
  ggpar(p, legend.title = "VAF")

  ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/MetaAML_vaf_distribution_median_corrected_raincloud.pdf", dpi = 300, width = 10, height = 2.5, units = "in")

  write.csv(threshold_list_all, "~/Desktop/MetaAML_results/Figure_3/MetaAML_vaf_distribution_median_corrected_raincloud.csv")


# analyze the clinical features that associate with high or low vaf per gene
sub1 = subset(final_data_matrix, mut_freq_gene > 50 & final_data_matrix$Subset == "de_novo" & final_data_matrix$Gene != "MLL")

# make sure that the FLT3 symbols are annotated well
for(i in 1:nrow(sub1)){
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] %in% c("ITD", "INDEL")){
    sub1$Gene[i] <- "FLT3-ITD"
  }
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] %in% c("SNV", "Deletion", "other")){
    sub1$Gene[i] <- "FLT3-TKD"
  }
}

genes = data.frame(unique(sub1$Gene))

# create directories for binary comparisions
# binary
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_Age")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_WBC")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_Platelet")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_LDH")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_Hemoglobin")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_PB_blast_percent")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/VAF_BM_blast_percent")

# summary plots for each variable ####
sub1$Gene <- as.factor(sub1$Gene)

variables = c("Age", "Platelet", "Hemoglobin", "WBC", "LDH", "BM_blast_percent", "PB_blast_percent")

colors = c("#bdbdbd",  "#fb6a4a", "#fe9929", "#d4b9da", "#238443", "#e6ab02", "#8dd3c7")

# continuous and discrete correlations ####
variables = data.frame("Age", "Platelet", "Hemoglobin", "WBC", "LDH", "BM_blast_percent", "PB_blast_percent")

variable_vaf_t_list = list()
y = 1

variable_vaf_list = list()
z = 1

for(i in 1:ncol(variables)){
  sub = sub1
  variable = as.character(variables[1,i])
  # var_feature = as.character(variables[1,i])
  
  print(variable)
  
  if(variable == "Platelet"){
    y_label <- "Platelet (1e3/μL)"
    y_max <- 375
    n_label <- 300
    point_color = "#fdd49e"
  } else if(variable == "Hemoglobin"){
    y_label <- "Hemoglobin (g/dL)"
    y_max <- 15
    n_label <- 12
    point_color = "#fe9929"
  } else  if(variable == "WBC"){
    y_label <- "WBC"
    y_max <- 300
    n_label = 240
    point_color = "#d4b9da"
  } else if(variable == "LDH"){
    y_label <- "LDH"
    y_max <- 1750
    n_label = 1400
    point_color = "#238443"
  } else if(variable == "Age"){
    y_label <- "Age"
    y_text = 90
    y_max = 100
    n_label = 80
    point_color = "#8073ac"
  } else if(variable == "BM_blast_percent"){
    y_label <- "BM Blast %"
    n_label = 80
    y_max = 109
    point_color = "#e6ab02"
  } else if(variable == "PB_blast_percent"){
    y_label <- "PB Blast %"
    n_label = 80
    y_max = 109
    point_color = "#8dd3c7"
  }
  
  sub$VAF <- as.numeric(sub$VAF)
  sub$VAF_CN_corrected <- as.numeric(sub$VAF_CN_corrected)
  sub$Age <- as.numeric(sub$Age)
  sub$Hemoglobin <- as.numeric(sub$Hemoglobin)
  sub$Platelet <- as.numeric(sub$Platelet)
  sub$WBC <- as.numeric(sub$WBC)
  sub$LDH <- as.numeric(sub$LDH)
  sub$BM_blast_percent <- as.numeric(sub$BM_blast_percent)
  sub$PB_blast_percent <- as.numeric(sub$PB_blast_percent)

  sub$median_vaf = as.factor(ifelse(sub$VAF_CN_corrected > median(sub$VAF_CN_corrected, na.rm = T), "Over", "Under"))
  sub = sub[!is.na(sub$median_vaf),]

  sub$median_vaf <- factor(sub$median_vaf , levels=c("Under", "Over"))
  
  cnum = which( colnames(sub)==variable)
  colnames(sub)[cnum] = "var_ft"
  
  p = ggplot(sub, aes(x=median_vaf, y=var_ft)) +
        geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
        geom_jitter(aes(fill = median_vaf), color = "black", shape = 21, position=position_jitter(0.2), size = 3) +
        scale_fill_manual(values = c("Over"=point_color, "Under" = "lightgrey")) +
        stat_compare_means(label.x = 1.1, label.y = y_max) +
        theme_cowplot(font_size = 20) +
        # labs(title = paste(gene)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab(label= variable) +
        xlab(label = NULL) +
        ylim(0,y_max) +
        theme(legend.position="right")  + labs(fill = "Median VAF") + theme(
          axis.text.x = element_blank(),
          axis.ticks.x.bottom  = element_blank())
  #
  p + facet_wrap(. ~ Gene, ncol = 6) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=1.5, linetype="solid"))
  #
  ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/", "VAF_", variable,"/",variable,"_discrete.png", sep = ""), dpi = 300, width = 15, height = 15, units = "in")
  
  
  # continuous
  # ggscatter(sub, "VAF", y = "Variable",
  #           color = "black", fill = point_color, size = 5, shape = 21,
  #           font.label = list(color = "black", size = 9, vjust = 0.5),
  #           # title = paste(gene),
  #           xlab = "VAF",
  #           ylab = y_label) +
  #   # xlim(0,1) +
  #   # ylim(0,y_max)+
  #   geom_smooth(method = "lm", color = "black", alpha = .5) +
  #   theme(plot.title = element_text(hjust = 0.5, size = 18), axis.title = element_text(size = 18)) +
  #   stat_cor(
  #     aes(label = paste(..rr.label..)),
  #     label.x.npc = 0.65,
  #     label.y.npc = 1
  #   ) +
  #   stat_cor(
  #     aes(label = paste(..p.label..)),
  #     label.x.npc = 0.65,
  #     label.y.npc = 0.9
  #   ) +
  #   # geom_text()
  #   annotate("text", label = paste("n =", n), x = 0.75, y = n_label, size = 4, colour = "black")
  
  # p + facet_wrap(. ~ Gene, ncol = 6) +
  #   theme(strip.background = element_rect(colour="black", fill="white", 
  #                                         size=1.5, linetype="solid"))
  
  # ggsave(filename = paste("~/Desktop/MetaAML_results/Data/Figures/", "VAF_", variable,"/",variable,".png", sep = ""), dpi = 300, width = 15, height = 15, units = "in") 
  
  
  
  for(j in 1:nrow(genes)){
    print(j)
    
    variable = as.character(variables[1,i])
    
    gene = as.character(genes$unique.sub1.Gene.[j])
    
    sub_mut = sub
    
    sub_mut <- subset(sub, sub$Gene == gene)
    
    sub_mut$VAF <- as.numeric(sub_mut$VAF)
    sub_mut$VAF_CN_corrected <- as.numeric(sub_mut$VAF_CN_corrected)
    sub_mut$Age <- as.numeric(sub_mut$Age)
    sub_mut$Hemoglobin <- as.numeric(sub_mut$Hemoglobin)
    sub_mut$Platelet <- as.numeric(sub_mut$Platelet)
    sub_mut$WBC <- as.numeric(sub_mut$WBC)
    sub_mut$LDH <- as.numeric(sub_mut$LDH)
    sub_mut$BM_blast_percent <- as.numeric(sub_mut$BM_blast_percent)
    sub_mut$PB_blast_percent <- as.numeric(sub_mut$PB_blast_percent)
  
    n <- as.numeric(nrow(sub_mut))
    print(n)
    
    names(sub_mut)[names(sub_mut) == "var_ft"] <- "Variable"
  
    # t-test on distribution of clinical variable based on VAF threshold 
    # find if the sample has a VAF greater or less than the median for that mutation  
    sub_mut$median_vaf = as.factor(ifelse(sub_mut$VAF_CN_corrected > median(sub_mut$VAF_CN_corrected, na.rm = T), "Over", "Under"))
    sub_mut = sub_mut[!is.na(sub_mut$median_vaf),]
    
    # calculate a p-value and effect size for the difference in clinical feaures based on VAF differences
    p_val =  wilcox.test(sub_mut$Variable ~ sub_mut$median_vaf, alternative = "two.sided")$p.value
    
    effect_size = cohens_d(data = sub_mut, Variable ~ median_vaf)$effsize
    n_n1 = as.numeric(length(which(sub_mut$median_vaf == "Over")))
    n_n2 = as.numeric(length(which(sub_mut$median_vaf == "Under")))
    effect_size_ci = psych::cohen.d.ci(d = effect_size, n = n, n1 = n_n1, n2 = n_n2)
    
    variable_vaf_t <- data.frame(matrix(NA, nrow = 1, ncol = 6))
    names(variable_vaf_t) <- c("Variable", "Mutated_Gene", "effect_size", "CI_lower", "CI_upper", "p_value")
    
    variable_vaf_t[1,1] <- paste(variable)
    variable_vaf_t[1,2] <- gene
    variable_vaf_t[1,3] <- effect_size
    variable_vaf_t[1,4] <- effect_size_ci[1]
    variable_vaf_t[1,5] <- effect_size_ci[3]
    variable_vaf_t[1,6] <- p_val
    
    # Add each list in the loop to a list of lists
    variable_vaf_t_list[[y]] = variable_vaf_t 
    
    y = y + 1
    
    # put the results of the t-test in a dataframe
    sub_mut$median_vaf <- factor(sub_mut$median_vaf , levels=c("Under", "Over"))
    
    ylab_max = max(sub_mut$Variable, na.rm = TRUE)
    
    ggplot(sub_mut, aes(x=median_vaf, y=Variable)) + 
      geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
      geom_jitter(aes(fill = median_vaf), color = "black", shape = 21, position=position_jitter(0.2), size = 3) +
      scale_fill_manual(values = c("Over"=point_color, "Under" = "lightgrey")) +
      stat_compare_means(label.x = 1.1, label.y = ylab_max) +
      theme_cowplot(font_size = 20) +
      labs(title = paste(gene)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylab(label= variable) +
      xlab(label = NULL) +
      ylim(0,ylab_max) +
      theme(legend.position="right")  + labs(fill = "Median VAF") + theme(
        axis.text.x = element_blank(),
        axis.ticks.x.bottom  = element_blank())
    
    ggsave(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/", "VAF_", variable,"/",variable, "_vs_", gene,"_VAF.png", sep = ""), dpi = 300, width = 5, height = 5, units = "in") 
    
  }   
}

# compile dataframes for discrete and continuous analyses
variable_vaf_comparison_discrete <- do.call(rbind, variable_vaf_t_list)

# correct for multiple hypotheses by variable
var2_adj = list()
a = 1

for(i in 1:ncol(variables)){
  print(i)
  variable = as.character(variables[1,i])
  
  var2= subset(variable_vaf_comparison_discrete, variable_vaf_comparison_discrete$Variable == variable)
  
  var2$q_val = p.adjust(var2$p_value)
  
  var2_adj[[a]] = var2 
  
  a = a + 1
  
}
var2_adj_list <- do.call(rbind, var2_adj)
# variable_vaf_comparison$z_score <- qnorm(variable_vaf_comparison$p_value)

# var2_adj_list$sig = as.factor(ifelse(var2_adj_list$q_val < 0.2, "q < 0.2", "q > 0.2"))

for(i in 1:nrow(var2_adj_list)){
  if(var2_adj_list$q_val[i] > 0.2){
    var2_adj_list$sig[i] = "q > 0.2"
  }
  if(var2_adj_list$q_val[i] < 0.2 & var2_adj_list$q_val[i] > 0.1){
    var2_adj_list$sig[i] = "q < 0.2"
  }
  if(var2_adj_list$q_val[i] < 0.1 & var2_adj_list$q_val[i] > 0.05){
    var2_adj_list$sig[i] = "q < 0.1"
  }
  if(var2_adj_list$q_val[i] < 0.05){
    var2_adj_list$sig[i] = "q < 0.05"
  }
}


var2_adj_list$Variable = gsub("BM_blast_percent", "BM blast %", var2_adj_list$Variable)
var2_adj_list$Variable = gsub("PB_blast_percent", "PB blast %", var2_adj_list$Variable)

# order the factors
var2_adj_list$Variable = factor(var2_adj_list$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM blast %', 'PB blast %', 'Age'))

# plot the discrete results
p = ggplot(var2_adj_list, aes(x = Mutated_Gene, y = effect_size, label = p_value)) +
  geom_hline(yintercept=0, linetype = "dashed", color = "black") +
  theme_cowplot() +
  geom_pointrange(size = 0.75, stat = "identity", alpha = 0.75,
                  aes(x = Mutated_Gene, ymin = CI_lower, ymax = CI_upper, y = effect_size, color = sig)) +
  scale_color_manual(values = c("q < 0.2" = "#FDE725FF", "q < 0.1" = "#1F968BFF", "q < 0.05" = "#482677FF", "q > 0.2" = "grey")) +
  ylab("Effect Size (high VAF vs. low VAF)")+
  xlab("")+
  theme(legend.position = "right", legend.title = element_blank(),
        axis.title.y=element_blank()) +
  coord_flip() +
  xlim(rev(levels(factor(var2_adj_list$Mutated_Gene))))

p + facet_grid(. ~ Variable)

p + facet_grid(. ~ Variable) +
  theme(strip.background = element_rect(colour="black", fill= "white",
                                        size=1.5, linetype="solid"))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/VAF_clinical_features_discrete.pdf", dpi = 300, width = 15, height = 6, units = "in") 

write.csv(var2_adj_list, "~/Desktop/MetaAML_results/Figure_3/VAF_discrete_clinical_features.csv")



# optimal vaf cutoff ####
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_vaf_thresholds")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_vaf_thresholds/logrank_plots")

load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")

final_data_matrix_2_sub = subset(final_data_matrix_2, final_data_matrix_2$mut_freq_gene >= 50 & final_data_matrix_2$Gene != "MLL" & final_data_matrix_2$Subset == "de_novo")
final_data_matrix_2_sub$Time_to_OS <- (final_data_matrix_2_sub$Time_to_OS/365)

final_data_matrix_2_sub$Gene = as.character(final_data_matrix_2_sub$Gene)

for(i in 1:nrow(final_data_matrix_2_sub)){
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] %in% c("SNV", "Deletion", "other")){
    final_data_matrix_2_sub$Gene[i] = "FLT3-TKD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] %in% c("ITD", "INDEL")){
    final_data_matrix_2_sub$Gene[i] = "FLT3-ITD"
  }
}

genes = data.frame(unique(final_data_matrix_2_sub$Gene))
genes = genes %>%
  arrange(genes,unique.final_data_matrix_2_sub.Gene.)
# now define the optimal VAF cutoff for each gene
results_list = list()
n = 1

plot_list = list()

for(i in 1:nrow(genes)){
  print(i)
  gene <- as.character(genes[i,1])
  
  final_data_matrix_2_sub2 <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$Gene == gene)
  final_data_matrix_2_sub2 = na.omit(distinct(final_data_matrix_2_sub2, Sample, Gene, VAF, VAF_CN_corrected, Time_to_OS, Censor))
  
  final_data_matrix_2_sub2 <- final_data_matrix_2_sub2[order(final_data_matrix_2_sub2$Sample, -final_data_matrix_2_sub2$VAF_CN_corrected),]
  final_data_matrix_2_sub2 = final_data_matrix_2_sub2[!duplicated(final_data_matrix_2_sub2$Sample),]
  
    if( nrow(final_data_matrix_2_sub2) >= 15){
      
    # find the cutoff
    final_data_matrix_2_sub2$vaf_threshold <- NA
    
    final_data_matrix_2_sub2$Censor = as.numeric(final_data_matrix_2_sub2$Censor)
    
    mstat <- maxstat.test(Surv(Time_to_OS, Censor) ~ VAF_CN_corrected, data=final_data_matrix_2_sub2, 
                          smethod="LogRank", pmethod="exactGauss", 
                          abseps=0.01)
    
    png(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_vaf_thresholds/logrank_plots/", gene, "_vaf_logrank_plot.png", sep = ""),res = 300, width = 5, height = 5, units = "in")
    plot(mstat)
    dev.off() 
    
    threshold = round(as.numeric(mstat$estimate), 2)
    
    for(i in 1:nrow(final_data_matrix_2_sub2)){
      if(!is.na(final_data_matrix_2_sub2$VAF_CN_corrected[i])){
        if(final_data_matrix_2_sub2$VAF_CN_corrected[i] >= threshold){
          final_data_matrix_2_sub2$vaf_threshold[i] = "over"
        }
        if(final_data_matrix_2_sub2$VAF_CN_corrected[i] < threshold){
          final_data_matrix_2_sub2$vaf_threshold[i] = "under"
        }
      }
    }
    
    
    n_over = length(which(final_data_matrix_2_sub2$vaf_threshold == "over"))
    n_under  = length(which(final_data_matrix_2_sub2$vaf_threshold == "under"))
    
    if(n_distinct(final_data_matrix_2_sub2$vaf_threshold) > 1){
      
      # summarize results from a Cox model
      final_data_matrix_2_sub2$vaf_threshold = ifelse(final_data_matrix_2_sub2$vaf_threshold == "over", 1,0)
      
      OS_fit <- survfit(Surv(Time_to_OS, Censor) ~ 1, data=final_data_matrix_2_sub2)
      OS_trt_fit <- survfit(Surv(Time_to_OS, Censor) ~ vaf_threshold, data=final_data_matrix_2_sub2, conf.type = "log-log")
      
      # run the Cox model
      model <- coxph( Surv(Time_to_OS, Censor) ~ vaf_threshold,
                      data = final_data_matrix_2_sub2 )
      
      # extract the informative data from the survival model
      array_dat = summary(model)$conf.int[1:4]
      array_dat[5] = gene
      
      # extract the log-rank p-value for the individual comparisons
      array_dat[6] = summary(model)$sctest[3]
      array_dat = array_dat[-2]
      
      n_pts = ifelse(array_dat[1] < 1, n_under, n_over)
      
      forest_plot_data <- data.frame("Gene" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = round(as.numeric(array_dat[2]),2), "Upper_CI" = round(as.numeric(array_dat[3]),2), "log_rank_p" = array_dat[5], "VAF_threshold" = threshold, "n_pts" = n_pts)
      
      results_list[[n]] <- forest_plot_data
      n=n+1   
      
      p = round(summary(model)$sctest[3],3)
      
      if(p <= 0.05){
        print(gene)
        
        p = ifelse(p < 0.001, paste0("p < 0.001"), paste("p =", p))
        
        hr = paste("HR = ", round(as.numeric(forest_plot_data$HR), 2), " (", forest_plot_data$Lower_CI, "-", forest_plot_data$Upper_CI, ")", sep = "")
        
        p_hr = paste(p, "\n", hr, sep = "")
        
        cohorts <- c("under", "over")
        
        title <- paste(threshold, "VAF\nthreshold",sep = " ")
        
        # plots the survival
        plot_list[[i]] = ggsurvplot(OS_trt_fit,
                                data = final_data_matrix_2_sub2,
                                log = (OS_trt_fit),
                                log.rank.weights = c("survdiff"),
                                pval = paste0(p_hr),
                                test.for.trend = F,
                                pval.method.size = 3,
                                pval.coord = c(3, 0.95),
                                conf.int = F,
                                censor = T,
                                surv.median.line = "none",
                                risk.table = T,
                                risk.table.title = c("Number at risk"),
                                risk.table.fontsize = 4,
                                risk.table.height = .3,
                                risk.table.y.text = T,
                                break.time.by = 5,
                                risk.table.pos = c("out"),
                                palette = c("#80796BFF", "#6A6599FF"),
                                xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(15, "plain", "#252525"),
                                ggtheme = theme_cowplot(),
                                pval.size = 4,
                                font.x = c(15, "plain", "#252525"),
                                font.y =  c(15, "plain", "#252525"),
                                font.legend = c(15, "plain"),
                                font.tickslab = c(15, "plain", "#252525"),
                                legend.labs = cohorts,
                                legend.title = paste(title),
                                legend = "right",
                                title = gene)
        print(plot_list[[i]])
        png(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_vaf_thresholds/", gene, ".png", sep = ""), res = 300, width = 5, height = 5, units = "in")
        #
        plot_list[[i]]
        print(plot_list[[i]])
        dev.off()
      }
    }
  }
}

temp_final = as.data.frame(do.call(rbind, results_list))
temp_final$q_value = p.adjust(temp_final$log_rank_p, method = "fdr")

write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_vaf_thresholds/optimal_vaf_threshold_for_survival_prediction.csv", row.names = F)

# summary plot of all optimal survival curves
plot_list = list.clean(plot_list)

plots = arrange_ggsurvplots(plot_list, print = FALSE,
                            ncol = 6, nrow = 2)
ggsave("~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_VAF_survival_grid.pdf", plots, width = 25, height = 10)

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
temp_final_hr$q_text = NA
for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$fdr[i] < 0.2){
    temp_final_hr$q_text[i] = temp_final_hr$fdr[i]
  }
}

temp_final_hr$p_text = as.numeric(temp_final_hr$p_text)
temp_final_hr$q_text = as.numeric(temp_final_hr$q_text)

temp_final_hr$p_text = round(temp_final_hr$p_text, 3)
temp_final_hr$q_text = round(temp_final_hr$q_text, 2)

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] < 0.01){
    temp_final_hr$p_text[i] = "p < 0.01"
  }
  if(temp_final_hr$log_rank_p[i] >= 0.01 & temp_final_hr$log_rank_p[i] <= 0.05){
    temp_final_hr$p_text[i] = paste("p =", paste(temp_final_hr$p_text[i]))
  }
}

temp_final_hr$p_q_text = paste(temp_final_hr$p_text, "; q = ", temp_final_hr$q_text, sep = "")

for(i in 1:nrow(temp_final_hr)){
  if(temp_final_hr$log_rank_p[i] > 0.05){
    temp_final_hr$p_q_text[i] = ""
  }
  if(temp_final_hr$log_rank_p[i] > 0.05){
    temp_final_hr$p_q_text[i] = ""
  }
}

temp_final_hr$HR = as.numeric(temp_final_hr$HR)
temp_final_hr$Lower_CI = as.numeric(temp_final_hr$Lower_CI)
temp_final_hr$Upper_CI = as.numeric(temp_final_hr$Upper_CI)

temp_final_hr$Gene <- factor(temp_final_hr$Gene, levels = temp_final_hr$Gene[order(temp_final_hr$HR)])

ggplot(temp_final_hr, aes(x = reorder(Gene, -HR), y = HR, label = p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(Gene, Upper_CI), hjust = 0, nudge_y = 0.5) +
  geom_point(aes(size = n_pts, color = sig_color), shape = 19, alpha = .9) +
  geom_segment(data = temp_final_hr, aes(y = Lower_CI, yend = Upper_CI, x = Gene, xend = Gene, col= sig_color), size = 0.5) +
  scale_color_manual(values = c("0" = "grey", "1" = "#1b7837", "2" = "#762a83"))+
  ylab("Hazard Ratio\n(high VAF vs. low VAF)")+
  theme_cowplot() +
  scale_size_area(max_size = 5,breaks=c(25,50,100,200,300)) +
  theme(legend.position = c(0.6, .2),
        axis.title.y=element_blank()) +
  coord_flip(ylim = c(0,22)) +
  guides(color = FALSE,
         size =  guide_legend(override.aes=list(colour="lightgrey"), title = "n. patients over\nVAF threshold"))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/gene_vaf_optimal_vaf_thresholds_hr_forest_plot_de_novo.pdf", dpi = 300, width = 5, height = 6, units = "in")


# HR differences ####
gene_hrs = read.csv("~/Desktop/MetaAML_results/Figure_2/Supplimental/gene_hazard_ratio.csv")
gene_vaf_hrs = read.csv("~/Desktop/MetaAML_results/Figure_3/Supplimental/optimal_vaf_thresholds/optimal_vaf_threshold_for_survival_prediction.csv")

gene_hrs = gene_hrs %>% 
  select(gene, HR, lower_95, upper_95, q_value)

names(gene_hrs) = c("Gene", "HR", "gene_lower_95", "gene_upper_95", "gene_q_value")

gene_vaf_hrs = gene_vaf_hrs %>% 
  select(Gene, HR, Lower_CI, Upper_CI, q_value, VAF_threshold)

names(gene_vaf_hrs) = c("Gene", "VAF_HR", "VAF_lower_95", "VAF_upper_95", "VAF_q_value", "VAF_threshold")

hr_comparision = left_join(gene_vaf_hrs, gene_hrs, by = "Gene")

hr_comparision$sig = ifelse(hr_comparision$VAF_q_value < 0.1, "q < 0.1", "q > 0.1")

hr_comparision$ratio = hr_comparision$VAF_HR/hr_comparision$HR
hr_comparision$label = NA

for(i in 1:nrow(hr_comparision)){
  if(hr_comparision$VAF_q_value[i] <= 0.1 | hr_comparision$ratio[i] > 2 | hr_comparision$ratio[i] < 0.5){
    hr_comparision$label[i] = hr_comparision$Gene[i]
  }
}


p =  ggplot(hr_comparision, aes(x=HR, y = VAF_HR, color = sig)) +
        theme_cowplot() +
  geom_segment(aes(x = 0, xend =4, y = 1, yend = 1),
               colour = "lightgrey",lty = "dashed") +
        geom_point(aes(size = VAF_threshold), alpha = .9) +
  geom_point(aes(size = VAF_threshold), shape = 1, colour = "black") +
        scale_colour_manual(values = c("q < 0.1"= "#b2182b","q > 0.1"="grey")) + 
        geom_label_repel(aes(label=label),size = 3, force = 25, show_guide = F, max.overlaps = 15) +
        coord_cartesian(ylim=c(0,5.5), xlim=c(0,4)) +
        ylab(label= "Hazard Ratio\n(high VAF vs. low VAF)") +
        xlab(label= "Hazard Ratio\n(Mutated vs. WT)") +
        theme(legend.position=c(0.65,0.75)) 
  
g = guide_legend(override.aes=list(colour="grey"), "VAF threshold")
h = guide_legend("VAF HR")
p  + guides(color = h, size = g)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/HR_vs_VAF_HR.pdf", dpi = 300, width = 4.5, height = 6, units = "in")
write_csv(hr_comparision, "~/Desktop/MetaAML_results/Figure_3/HR_vs_VAF_HR.csv")


# Effect size differences ####
gene_effect_size = read.csv("~/Desktop/MetaAML_results/Figure_2/Supplimental/gene_clinical_features.csv")
gene_vaf_effect_size = read.csv("~/Desktop/MetaAML_results/Figure_3/VAF_discrete_clinical_features.csv")

gene_effect_size = gene_effect_size %>% 
  select(Variable, Gene, effect_size, n_mut, q_val) %>%
  subset(Variable %in% c("WBC", "PB_blast_percent"))
gene_effect_size$Variable = gsub("PB_blast_percent", "PB blast %", gene_effect_size$Variable)

names(gene_effect_size) = c("Variable", "Gene", "effect_size", "n_mut", "q_value")

gene_vaf_effect_size = gene_vaf_effect_size %>% 
  select(Variable, Mutated_Gene, effect_size, q_val) %>%
  subset(Variable %in% c("WBC", "PB blast %"))

names(gene_vaf_effect_size) = c("Variable", "Gene", "effect_size_vaf", "q_value_vaf")

es_comparision = inner_join(gene_vaf_effect_size, gene_effect_size, by = c("Variable", "Gene"))

es_comparision$sig_color = "q > 0.1"

for(i in 1:nrow(es_comparision)){
  if(es_comparision$q_value_vaf[i] < 0.1 & es_comparision$Variable[i] == "WBC"){
    es_comparision$sig_color[i] = "WBC q < 0.1"
  }
  if(es_comparision$q_value_vaf[i] < 0.1 & es_comparision$Variable[i] == "PB blast %"){
    es_comparision$sig_color[i] = "PB q < 0.1"
  }
}

es_comparision$ratio = es_comparision$effect_size_vaf/es_comparision$effect_size
es_comparision$label = NA

for(i in 1:nrow(es_comparision)){
  if(es_comparision$q_value_vaf[i] < 0.1){
    es_comparision$label[i] = es_comparision$Gene[i]
  }
}

es_comparision$Variable = factor(es_comparision$Variable, levels = c("WBC", "PB blast %"))

p = ggplot(es_comparision, aes(x=effect_size, y = effect_size_vaf, color = factor(sig_color))) +
  theme_cowplot() +
  geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
  geom_vline(xintercept = 0,  linetype = "dashed", color = "lightgrey") +
  geom_point(aes(size = n_mut), alpha = 1) +
  geom_point(aes(size = n_mut), shape = 1, colour = "black") +
  scale_colour_manual(values = c("WBC q < 0.1"= "#d4b9da","PB q < 0.1"="#8dd3c7", "q > 0.1" = "grey")) + 
  geom_label_repel(aes(label=label),size = 3, force = 50, show_guide = F, max.overlaps = 10) +
  ylab(label= "Effect Size\n(high VAF vs. low VAF)") +
  xlab(label= "Effect Size\n(Mutated vs. WT)") +
  theme(legend.position="right") 

g = guide_legend(override.aes=list(colour="grey"), "n. patients")
h = guide_legend("VAF effect size")

p + facet_wrap(. ~ Variable, ncol = 1, nrow = 2) + guides(size = g, color = FALSE) +
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid")) + guides(color = h, size = g) 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/WBC_ES_vs_VAF_ES.pdf", dpi = 300, width = 5, height = 6, units = "in")
write_csv(es_comparision, "~/Desktop/MetaAML_results/Figure_3/WBC_ES_vs_VAF_ES.csv")


# Supplimental ####

# static vaf cutoff ####
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/Median_VAF")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental/30_VAF")

load("~/Desktop/MetaAML_results/final_data_matrix_2.RData")

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
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "Deletion"){
    final_data_matrix_2_sub$Gene[i] = "FLT3-TKD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "INDEL"){
    final_data_matrix_2_sub$Gene[i] = "FLT3-ITD"
  }
  if(final_data_matrix_2_sub$Gene[i] == "FLT3" & final_data_matrix_2_sub$variant_type[i] == "other"){
    final_data_matrix_2_sub$Gene[i] = "FLT3-TKD"
  }
}

n=n_distinct(final_data_matrix_2_sub$Gene)
genes=data.frame(unique(final_data_matrix_2_sub$Gene))


results_list = list()
results_list_median = list()

n=1
for(i in 1:nrow(genes)){
  gene=genes[i,1]
  mut_pts=subset(final_data_matrix_2_sub, Gene == genes[i,1])
  mut_pts = distinct(mut_pts, Sample, Gene, VAF, VAF_CN_corrected, Time_to_OS, Censor)
  mut_pts$hr_stratifier_vaf = 0
  mut_pts$hr_stratifier_vaf_text = "Low"
  mut_pts$hr_stratifier_vaf_median = 0
  mut_pts$hr_stratifier_vaf_median_text = "Below"
  
  mut_pts <- mut_pts[order(mut_pts$Sample, -mut_pts$VAF_CN_corrected),]
  mut_pts= mut_pts[!duplicated(mut_pts$Sample),]
  
  for(j in 1:nrow(mut_pts)){
    if(!is.na(mut_pts$VAF_CN_corrected[j])){
      # try splitting by median and 0.3
      if(mut_pts$VAF_CN_corrected[j] >= 30){
        mut_pts$hr_stratifier_vaf[j] = 1
        mut_pts$hr_stratifier_vaf_text[j] = "High"
      } 
      if(mut_pts$VAF_CN_corrected[j] >= mean(mut_pts$VAF_CN_corrected, na.rm = T)){
        mut_pts$hr_stratifier_vaf_median[j] = 1
        mut_pts$hr_stratifier_vaf_median_text[j] = "Above"
      } 
    }
  }
  
  # run the Cox model
  mut_pts$Time_to_OS=as.numeric(mut_pts$Time_to_OS)
  
  mut_pts$Censor = as.numeric(mut_pts$Censor)
  
  # analyze the dynamic median VAF survival results
  model_median <- coxph( Surv(Time_to_OS, Censor) ~ hr_stratifier_vaf_median,
                         data = mut_pts)
  
  # extract the informative data from the survival model
  array_dat = summary(model_median)$conf.int[1:4]
  array_dat[5] = genes[i,1]
  
  # extract the log-rank p-value for the individual comparisons
  array_dat[6] = summary(model_median)$sctest[3]
  array_dat = array_dat[-2]
  
  forest_plot_data_median <- data.frame("Gene" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = array_dat[2], "Upper_CI" = array_dat[3], "log_rank_p" = array_dat[5])
  
  results_list_median[[n]] <- forest_plot_data_median
  
  p1=forest_plot_data_median$log_rank_p
  
  # plot if significant
  if(p1 <= 0.05){
    # plots the survival
    
    OS_fit <- survfit(Surv(Time_to_OS, Censor) ~ 1, data=mut_pts)
    OS_trt_fit <- survfit(Surv(Time_to_OS, Censor) ~ hr_stratifier_vaf_text, data=mut_pts, conf.type = "log-log")
    
    surv_plot <- ggsurvplot(OS_trt_fit,
                            data = mut_pts,
                            log = (OS_trt_fit),
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
                            risk.table.height = .4,
                            risk.table.y.text = T,
                            break.time.by = 5,
                            risk.table.pos = c("out"),
                            palette = c("Above" = "#6A6599FF", "Below" = "#80796BFF"),
                            xlab = "Years",
                            ylim = c(0, 1.0),
                            ylab =  "Survival Probability",
                            font.main = c(15, "plain", "#252525"),
                            pval.size = 4,
                            font.x = c(12, "plain", "#252525"),
                            font.y =  c(12, "plain", "#252525"),
                            font.legend = c(12, "plain"),
                            font.tickslab = c(12, "plain", "#252525"),
                            legend.labs = c("0" = "Above", "1" = "Below"),
                            legend.title = paste("Median VAF"),
                            legend = "right",
                            title = gene,
                            ggtheme = theme_cowplot()
                            # ggtheme = theme(plot.title = element_text(hjust = 0.5))
    )
    
    print(surv_plot)
    png(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/Median_VAF/",gene,"_survival_by_median_VAF.png", sep = ""), res = 300, width = 4, height = 4, units = "in")
    
    surv_plot
    print(surv_plot)
    dev.off()
  }
  
  
  
  ## analyze the 0.3 threshold survival results
  model <- coxph( Surv(Time_to_OS, Censor) ~ hr_stratifier_vaf,
                  data = mut_pts)
  
  # extract the informative data from the survival model
  array_dat = summary(model)$conf.int[1:4]
  array_dat[5] = genes[i,1]
  
  # extract the log-rank p-value for the individual comparisons
  array_dat[6] = summary(model)$sctest[3]
  array_dat = array_dat[-2]
  
  forest_plot_data <- data.frame("Gene" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = array_dat[2], "Upper_CI" = array_dat[3], "log_rank_p" = array_dat[5])
  
  results_list[[n]] <- forest_plot_data
  n=n+1   
  
  p2=forest_plot_data$log_rank_p
  
  # plot if significant
  if(p2 <= 0.05){
    # plots the survival
    
    OS_fit <- survfit(Surv(Time_to_OS, Censor) ~ 1, data=mut_pts)
    OS_trt_fit <- survfit(Surv(Time_to_OS, Censor) ~ hr_stratifier_vaf_text, data=mut_pts, conf.type = "log-log")
    
    surv_plot <- ggsurvplot(OS_trt_fit,
                            data = mut_pts,
                            log = (OS_trt_fit),
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
                            risk.table.height = .4,
                            risk.table.y.text = T,
                            break.time.by = 5,
                            risk.table.pos = c("out"),
                            palette = c("High" = "#6A6599FF", "Low" = "#80796BFF"),
                            xlab = "Years",
                            ylim = c(0, 1.0),
                            ylab =  "Survival Probability",
                            font.main = c(15, "plain", "#252525"),
                            pval.size = 4,
                            font.x = c(12, "plain", "#252525"),
                            font.y =  c(12, "plain", "#252525"),
                            font.legend = c(12, "plain"),
                            font.tickslab = c(12, "plain", "#252525"),
                            legend.labs = c("0" = "High", "1" = "Low"),
                            legend.title = paste("VAF"),
                            legend = "right",
                            title = gene,
                            ggtheme = theme_cowplot()
                            # ggtheme = theme(plot.title = element_text(hjust = 0.5))
    )
    
    print(surv_plot)
    png(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/30_VAF/",gene,"_survival_by_30_VAF.png", sep = ""), res = 300, width = 4, height = 4, units = "in")
    
    surv_plot
    print(surv_plot)
    dev.off()
  }
}

temp_final = as.data.frame(do.call(rbind, results_list))
temp_final_median = as.data.frame(do.call(rbind, results_list_median))

# plot the results for 0.3 VAF cuttoff
temp_final$Gene <- factor(temp_final$Gene, levels = temp_final$Gene[order(temp_final$HR)])

temp_final$sig_color = 0

for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] < 0.05){
    temp_final$sig_color[i] =1
  }
}

temp_final$sig_color = as.factor(temp_final$sig_color)

temp_final$fdr = p.adjust(temp_final$log_rank_p, method = "fdr")

temp_final$sig_color = as.factor(temp_final$sig_color)


temp_final = subset(temp_final, temp_final$Upper_CI != "Inf")

temp_final$p_text = NA
temp_final$q_text = NA
for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] < 0.05){
    temp_final$p_text[i] = temp_final$log_rank_p[i]
    temp_final$q_text[i] = temp_final$fdr[i]
  }
}

temp_final$p_text = as.numeric(temp_final$p_text)
temp_final$q_text = as.numeric(temp_final$q_text)

temp_final$p_text = round(temp_final$p_text, 3)
temp_final$q_text = round(temp_final$q_text, 2)

temp_final$p_q_text = paste("p = ", temp_final$p_text, "; q = ", temp_final$q_text, sep ="")
temp_final$p_text = paste("p = ", temp_final$p_text)

for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] > 0.05){
    temp_final$p_q_text[i] = ""
  }
  if(temp_final$log_rank_p[i] > 0.05){
    temp_final$p_text[i] = ""
  }
}

temp_final$HR = as.numeric(temp_final$HR)
temp_final$Upper_CI = as.numeric(temp_final$Upper_CI)
temp_final$Lower_CI = as.numeric(temp_final$Lower_CI)

ggplot(temp_final, aes(x = reorder(Gene, -HR), y = HR, label = temp_final$p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=3, linetype="dashed", color = "lightgrey") +
  geom_text(aes(Gene, Upper_CI), hjust = 0, nudge_y = 0.1) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 19, 
                  # fill = "white",
                  aes(x = Gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "#762a83"))+
  ylab("Hazard Ratio")+
  theme_cowplot() +
  ylim(0,17) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip() 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/30_VAF/gene_vaf_discrete_hr_forest_plot_de_novo_30.pdf", dpi = 300, width = 5, height = 7.5, units = "in")

write.csv(temp_final, "~/Desktop/MetaAML_results/Figure_3/Supplimental/30_VAF/static_vaf_threshold_survival.csv")


# plot the results for median VAF cuttoff
temp_final_median$Gene <- factor(temp_final_median$Gene, levels = temp_final_median$Gene[order(temp_final_median$HR)])

temp_final_median$sig_color = 0

for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$log_rank_p[i] < 0.05){
    temp_final_median$sig_color[i] =1
  }
}

temp_final_median$sig_color = as.factor(temp_final_median$sig_color)

temp_final_median$fdr = p.adjust(temp_final_median$log_rank_p, method = "fdr")

temp_final_median$sig_color = as.factor(temp_final_median$sig_color)

temp_final_median = subset(temp_final_median, temp_final_median$Upper_CI != "Inf")

temp_final_median$categories <- reorder(temp_final_median$categories, temp_final_median$HR)

temp_final_median$p_text = NA
for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$log_rank_p[i] < 0.05){
    temp_final_median$p_text[i] = temp_final_median$log_rank_p[i]
  }
}
temp_final_median$q_text = NA
for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$fdr[i] < 0.2){
    temp_final_median$q_text[i] = temp_final_median$fdr[i]
  }
}

temp_final_median$p_text = as.numeric(temp_final_median$p_text)
temp_final_median$q_text = as.numeric(temp_final_median$q_text)

temp_final_median$p_text = round(temp_final_median$p_text, 3)
temp_final_median$q_text = round(temp_final_median$q_text, 2)

for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$log_rank_p[i] < 0.01){
    temp_final_median$p_text[i] = "p < 0.01"
  }
  if(temp_final_median$log_rank_p[i] >= 0.01 & temp_final_median$log_rank_p[i] <= 0.05){
    temp_final_median$p_text[i] = paste("p = ", temp_final_median$p_text[i], sep = "")
  }
}

for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$fdr[i] <= 0.2){
    temp_final_median$p_q_text = paste(temp_final_median$p_text, "; q = ", temp_final_median$q_text, sep = "")
  }
}

for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$log_rank_p[i] > 0.05){
    temp_final_median$p_q_text[i] = ""
  }
  if(temp_final_median$log_rank_p[i] > 0.05){
    temp_final_median$p_text[i] = ""
  }
  if(temp_final_median$fdr[i] > 0.2){
    temp_final_median$p_q_text[i] = ""
  }
}


temp_final_median$HR = as.numeric(temp_final_median$HR)
temp_final_median$Upper_CI = as.numeric(temp_final_median$Upper_CI)
temp_final_median$Lower_CI = as.numeric(temp_final_median$Lower_CI)

ggplot(temp_final_median, aes(x = reorder(Gene, -HR), y = HR, label = temp_final_median$p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=3, linetype="dashed", color = "lightgrey") +
  geom_text(aes(Gene, Upper_CI), hjust = 0, nudge_y = 0.1) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 19, 
                  # fill = "white",
                  aes(x = Gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "darkgrey", "1" = "#762a83"))+
  ylab("Hazard Ratio")+
  theme_cowplot() +
  # ylim(0,18)+
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/Median_VAF/gene_vaf_discrete_hr_forest_plot_de_novo_median.pdf", dpi = 300, width = 5, height = 7.5, units = "in")

write.csv(temp_final_median, "~/Desktop/MetaAML_results/Figure_3/Supplimental/Median_VAF/static_vaf_threshold_survival.csv")



