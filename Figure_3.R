# ========================================================================================================================================== #
# Figure_3.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 09/03/2020
# Description: This script will perform survival analyses based on VAF thresholds as seen in Figure 3 of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# libraries
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
if (!require('maxstat')) install.packages('maxstat'); library('maxstat')
if (!require('rstatix')) install.packages('rstatix'); library('rstatix')
if (!require('psych')) install.packages('psych'); library('psych')

library(maxstat)

dir.create("~/Desktop/MetaAML_results/Figure_3")
dir.create("~/Desktop/MetaAML_results/Figure_3/Supplimental")

# VAf distribution ####
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

vaf_sub = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
vaf_sub = subset(vaf_sub, vaf_sub$mut_freq_gene >= 50)
vaf_sub = subset(vaf_sub, vaf_sub$Gene != "MLL")

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(vaf_sub)){
  if(vaf_sub$Gene[i] == "FLT3" & vaf_sub$variant_type[i] == "ITD"){
    vaf_sub$Gene[i] <- "FLT3-ITD"
  }
  if(vaf_sub$Gene[i] == "FLT3" & vaf_sub$variant_type[i] == "SNV"){
    vaf_sub$Gene[i] <- "FLT3-TKD"
  }
  if(vaf_sub$Gene[i] == "FLT3" & vaf_sub$variant_type[i] == "Deletion"){
    vaf_sub$Gene[i] <- "FLT3-TKD"
  }
  if(vaf_sub$Gene[i] == "FLT3" & vaf_sub$variant_type[i] == "INDEL"){
    vaf_sub$Gene[i] <- "FLT3-ITD"
  }
}

vaf_sub = select(vaf_sub, Sample, Gene, VAF_male_x, variant_type)

vaf_sub$threshold = ifelse(vaf_sub$VAF_male_x >= .3 , "High", "Low")

vaf_sub = subset(vaf_sub, vaf_sub$threshold == "High" | vaf_sub$threshold == "Low")

vaf_sub$VAF_male_x=as.numeric(vaf_sub$VAF_male_x)
vaf_sub$Gene = as.character(vaf_sub$Gene)

vaf_sub$Gene <- with(vaf_sub, reorder(Gene, -VAF_male_x, median))


p = ggplot(vaf_sub, aes(x=Gene, y=VAF_male_x)) + 
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
  geom_jitter(aes(fill = threshold), color = "black", shape = 21, position=position_jitter(0.2), size = 1.5) +
  scale_fill_manual(values = c("#cb181d", "#3690c0")) +
  # scale_fill_manual(values = c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")) +
  geom_hline(yintercept = .3, color = "#b2182b", linetype="dashed") +
  theme_cowplot(font_size = 7.5) +
  labs(title = NULL) +
  ylab(label = "VAF") +
  xlab(label = NULL) +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7.5, hjust = 1),
        axis.text.y = element_text(size = 7.5))

ggpar(p, legend.title = "")
ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/MetaAML_vaf_distribution.pdf", dpi = 300, width = 10, height = 3, units = "in")


# now plot the VAF distributino per gene colored by whethere the mutation is above or below the mdeian VAF for that gene
genes = unique(vaf_sub$Gene)

vaf_sub$median_threshold = NA

threshold_list = list()
z = 1

for(i in 1:length(genes)){
  # print(i)
  # select the gene of interest
  sub = subset(vaf_sub, Gene == genes[i])
  
  # calculate the median vaf based on the X-corrected VAF
  med_vaf = as.numeric(median(sub$VAF_male_x))
  
  sub$median_threshold = ifelse(sub$VAF_male_x >= med_vaf, "High", "Low")
  
  threshold_list[[i]] = data.frame(sub)
  
  # z = z + 1
}

threshold_list_all = do.call(rbind, threshold_list)
rm(threshold_list)

threshold_list_all$Gene <- with(threshold_list_all, reorder(Gene, -VAF_male_x, median))

p = ggplot(threshold_list_all, aes(x=Gene, y=VAF_male_x)) + 
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
  geom_jitter(aes(fill = median_threshold), color = "black", shape = 21, position=position_jitter(0.2), size = 1.5) +
  scale_fill_manual(values = c("#cb181d", "#3690c0")) +
  # scale_fill_manual(values = c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")) +
  # geom_point(data=hline, aes(tt, v), shape=95, size=5) +
  # geom_hline(yintercept = hline, color = "#b2182b", linetype="dashed") +
  theme_cowplot(font_size = 7.5) +
  labs(title = NULL) +
  ylab(label = "VAF") +
  xlab(label = NULL) +
  theme(legend.position="right", legend.text = element_text(size = 7.5), axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7.5, hjust = 1),
        axis.text.y = element_text(size = 7.5))

ggpar(p, legend.title = "")
ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/MetaAML_vaf_distribution_median.pdf", dpi = 300, width = 10, height = 2.5, units = "in")


# alternate visualization using Rainplot
if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('jorvlan/raincloudplots')

library(raincloudplots)

packages <- c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown", 
              "readr", "caTools", "bitops")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

p=ggplot(threshold_list_all, aes(x = Gene, y = VAF_male_x, fill = median_threshold)) +
  geom_flat_violin(aes(fill = "darkgrey"),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = "black")+
  geom_boxplot(aes(x = Gene, y = VAF_male_x, fill = "darkgrey"),outlier.shape = NA, alpha = .5, width = .25, colour = "black")+
  geom_point(aes(x = as.numeric(Gene)-.15, y = VAF_male_x, colour = median_threshold),position = position_jitter(width = .025), size = 1.5, shape = 21, color = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  # scale_fill_brewer(palette = "Dark2")+
  scale_color_manual(values = c("darkgrey", "#cb181d", "#3690c0")) +
  scale_fill_manual(values = c( "darkgrey", "#cb181d", "#3690c0")) +
  theme_cowplot(font_size = 7.5) +
  ylim(0,1) +
  ggtitle("") +
  ylab("VAF") +
  xlab("") +
  theme(axis.text.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7.5, hjust = 1))
ggpar(p, legend.title = "")
ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/MetaAML_vaf_distribution_median_detailed.pdf", dpi = 300, width = 10, height = 2.5, units = "in")


ggplot(threshold_list_all, aes(x = Gene, y = VAF_male_x, fill = median_threshold)) +
  # geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = "grey")+
  geom_jitter(aes(fill = median_threshold), color = "black", shape = 21, position=position_jitter(0.2), size = 1.5) +
  geom_boxplot(aes(x = Gene, y = VAF_male_x),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_fill_manual(values = c("#cb181d", "#3690c0")) +
  theme_cowplot() +
  ggtitle("") +
  ylab("") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))




#### CNA correction ####
# Using the cytogenetic data, recalculate the VAF for mutations falling on focal or general CNAs

# Pappamanuel/Gerstung VAF correction is performed using code provided by the authors
# TCGA VAFs are corrected using a gene-level CNA supplimentary file
# for all other cohorts, manual correction is performed using reported cytogenetic/karyotypic data. 

#### Papaemmanuil/Gerstung ####
download.file("https://github.com/gerstung-lab/AML-multistage/blob/master/data/AMLSG_Clinical_Anon.RData?raw=true", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")

load("~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")

# 1.3.2.2 Mutation data
download.file("https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/AMLSG_Genetic.txt", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Genetic.txt")

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
                Treatment = data.frame(ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL_os=tplIndexOs), 
                Demographics = clinicalData[,c("AOD","gender")],
                Clinical = cbind(clinicalData[, c("Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH","HB","platelet",
                                                  "Splenomegaly")], as.integer(clinicalData$TypeAML)),
                MolRisk = as.integer(clinicalData$M_Risk))
#
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
dataRaw = dataFrame
dim(dataFrame)
rownames(dataFrame) <- clinicalData$PDID


download.file("https://github.com/gerstung-lab/AML-multistage/raw/master/data/AMLSG_Karyotypes.txt", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Karyotypes.txt")

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

# plot the difference in VAF distribution per gene after CNA correction
raw_vaf = select(papaemmanuil_corrected, Gene, VAF)
raw_vaf$group = "VAF"

CNA_vaf = select(papaemmanuil_corrected, Gene, VAF_CN_corrected)
names(CNA_vaf)[2] = "VAF"
CNA_vaf$group = "VAF + \nCytogenetics"

all = rbind(raw_vaf, CNA_vaf)

# select mutations
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1", "FLT3-ITD","NRAS", "KRAS", "FLT3-TKD", "PTPN11")

all = subset(all, all$Gene %in% mut_list)

all$group <- factor(all$group , levels=c("VAF", "VAF + \nCytogenetics"))

p <- ggboxplot(all, x = "Gene", y = "VAF",
               color = "group", palette = "jama",
               add = "jitter",
               shape = 20)
p + stat_compare_means(aes(group = group), label = "p.signif") + xlab(NULL) +
  ylim(0,105) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + theme(legend.title=element_blank())






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

raw_vaf = select(tcga_dat, Gene, VAF)
raw_vaf$group = "VAF"

CNA_vaf = select(tcga_dat, Gene, VAF_CN_corrected)
names(CNA_vaf)[2] = "VAF"
CNA_vaf$group = "VAF + \nCytogenetics"

all = rbind(raw_vaf, CNA_vaf)


# select mutations
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1", "FLT3-ITD","NRAS", "KRAS", "FLT3-TKD", "PTPN11")

all = subset(all, all$Gene %in% mut_list)

all$group <- factor(all$group , levels=c("VAF", "VAF + \nCytogenetics"))

p <- ggboxplot(all, x = "Gene", y = "VAF",
               color = "group", palette = "jama",
               add = "jitter",
               shape = 20)
p + stat_compare_means(aes(group = group), label = "p.signif") + xlab(NULL) +
  ylim(0,105) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + theme(legend.title=element_blank())




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

# cohort_aggrigate$Cytogenetics = tolower(cohort_aggrigate$Cytogenetics)

# subset mutations to those where karyotyping data suggests deletions or amplifications
cohort_aggrigate$Cytogenetics = gsub('\\s+', '', cohort_aggrigate$Cytogenetics)
cohort_aggrigate$VAF_CN_corrected = NA

for(i in 1:nrow(cohort_aggrigate)){
  
  # find whole chromosome gains or losses
  ch_gain = paste("\\+",cohort_aggrigate$Chromosome.scaffold.name[i], ",", sep = "")
  ch_loss = paste("\\-",cohort_aggrigate$Chromosome.scaffold.name[i], ",", sep = "")
  
  # correct for VAF based on copy number gains for that chromosome  
  if(grepl(ch_gain, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$VAF[i] > 35) {cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*0.75}
  if(grepl(ch_gain, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$VAF[i] <= 35) {cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*1.5}
  # correct for VAF based on copy number loss for that chromosome
  if(grepl(ch_loss, cohort_aggrigate$Cytogenetics[i])){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2}
  
  # find focal chromosome gains or losses
  locus_gain = paste("add",cohort_aggrigate$partial_annotation[i],"|","add",cohort_aggrigate$full_annotation[i], sep = "")
  locus_loss = paste("del",cohort_aggrigate$partial_annotation[i],"|","del",cohort_aggrigate$full_annotation[i], sep = "")
  
  # correct for VAF based on broad copy number gains at the gene locus  
  if(grepl(locus_gain, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$VAF[i] > 35){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*0.75}
  if(grepl(locus_gain, cohort_aggrigate$Cytogenetics[i]) & cohort_aggrigate$VAF[i] <= 35){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*1.5}
  # correct for VAF based on broad copy number loss at the gene locus
  if(grepl(locus_loss, cohort_aggrigate$Cytogenetics[i])){ cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2}

}
for(i in 1:nrow(cohort_aggrigate)){
  print(i)
  if(is.na(cohort_aggrigate$VAF_CN_corrected[i])) { cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]}
}

# bind all of the cohorts together into a final data frame
# remove the copy number column from some dataframes to make it compatible with the other data frames
tcga_dat$cn = NULL
cohort_aggrigate[,28:32] = NULL

# bind all cohorts together
final_data_matrix = rbind(papaemmanuil_corrected, tcga_dat, cohort_aggrigate)
final_data_matrix_2 = final_data_matrix
save(final_data_matrix,  file = "~/Desktop/MetaAML_results/final_data_matrix_2.RData")


# plot the difference in VAF distribution per gene after CNA correction
raw_vaf = select(final_data_matrix, Gene, VAF, Cohort)
raw_vaf$group = "VAF"

CNA_vaf = select(final_data_matrix, Gene, VAF_CN_corrected, Cohort)
names(CNA_vaf)[2] = "VAF"
CNA_vaf$group = "VAF + \nCytogenetics"

all = rbind(raw_vaf, CNA_vaf)

# select mutations 
mut_list = c("TP53", "DNMT3A", "SRSF2", "IDH2", "JAK2", "TET2", "U2AF1", "IDH1", "ETV6", "RUNX1", "BCOR", "ASXL1", "PHF6", "SF3B1", "GATA2", "CBL", "CEBPA", "WT1", "RAD21", "EZH2", "NF1", "KIT", "NPM1", "FLT3-ITD","NRAS", "KRAS", "FLT3-TKD", "PTPN11")

all = subset(all, all$Gene %in% mut_list)

all$group <- factor(all$group , levels=c("VAF", "VAF + \nCytogenetics"))

# Box plot facetted by "gene"
p <- ggpaired(all, x = "group", y = "VAF",
              color = "group", palette = "jama", 
              line.color = "lightgray", line.size = 0.4,
              facet.by = "Gene", short.panel.labs = T)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE) + xlab(NULL) + ylab("VAF") + theme(legend.title=element_blank())
ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/all_vaf_correction_comparison_facet.pdf", dpi = 300, width = 10, height = 12, units = "in")

# subset to only those with VAF data and plot distributions by cohort
all = subset(all, all$Cohort %in% c("Papaemmanuil", "TCGA", "Tyner", "Majeti", "Welch", "Greif", "Hirsch"))

p <- ggboxplot(all, x = "Gene", y = "VAF",
               color = "group", palette = "jama",
               add = "jitter",
               # shape = "Cohort",
               facet.by = "Cohort",
               shape = 20)
p + stat_compare_means(aes(group = group), label = "p.signif") + xlab(NULL) +
  ylim(0,105) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + theme(legend.title=element_blank())

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/all_vaf_correction_comparison_facet.pdf", dpi = 300, width = 20, height = 10, units = "in")

# plot total vaf distribution for recurrent mutations
p <- ggboxplot(all, x = "Gene", y = "VAF",
               color = "group", palette = "jama",
               add = "jitter",
               # shape = "Cohort",
               # facet.by = "Cohort",
               shape = 20)
p + stat_compare_means(aes(group = group), label = "p.signif") + xlab(NULL) +
  ylim(0,105) +
  theme(legend.position="right", axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) + theme(legend.title=element_blank())

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/all_vaf_correction_comparison.pdf", dpi = 300, width = 15, height = 5, units = "in")


# now plot the VAF distribution per gene colored by if the mutation is above or below the median VAF for that gene
genes = unique(final_data_matrix$Gene)

final_data_matrix$median_threshold = NA

threshold_list = list()
z = 1

for(i in 1:length(genes)){
  # select the gene of interest
  sub = subset(final_data_matrix, Gene == genes[i])
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

p = ggplot(threshold_list_all, aes(x=Gene, y=VAF_CN_corrected)) + 
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
  geom_jitter(aes(fill = median_threshold_NA), color = "black", shape = 21, position=position_jitter(0.2), size = 1.5) +
  scale_fill_manual(values = c("#cb181d", "#3690c0")) +
  theme_cowplot(font_size = 7.5) +
  labs(title = NULL) +
  ylab(label = "VAF") +
  xlab(label = NULL) +
  # theme(legend.position="right") +
  theme(legend.position="right", legend.text = element_text(size = 7.5), axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                                                    size = 7.5, hjust = 1),
        axis.text.y = element_text(size = 7.5))

ggpar(p, legend.title = "")
ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/MetaAML_vaf_distribution_median_corrected.pdf", dpi = 300, width = 10, height = 2.5, units = "in")



# analyze the clinical features that associate with high or low vaf per gene
sub1 = subset(final_data_matrix, mut_freq_gene > 50 & final_data_matrix$Subset == "de_novo" & final_data_matrix$Gene != "MLL")

# make sure that the FLT3 symbols are annotated well
for(i in 1:nrow(sub1)){
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] == "ITD"){
    sub1$Gene[i] <- "FLT3-ITD"
  }
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] == "SNV"){
    sub1$Gene[i] <- "FLT3-TKD"
  }
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] == "Deletion"){
    sub1$Gene[i] <- "FLT3-TKD"
  }
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] == "INDEL"){
    sub1$Gene[i] <- "FLT3-ITD"
  }
  if(sub1$Gene[i] == "FLT3" & sub1$variant_type[i] == "other"){
    sub1$Gene[i] <- "FLT3-TKD"
  }
}

genes = data.frame(unique(sub1$Gene))

# create directories for binary comparisions
# binary
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Age/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_WBC/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Platelet/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_LDH/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_Hemoglobin/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_PB_blast_percent/discrete")
dir.create("~/Desktop/MetaAML_results/Data/Figures/VAF_BM_blast_percent/discrete")

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
        ylab(label= var_feature) +
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
  ggsave(filename = paste("~/Desktop/MetaAML_results/Data/Figures/", "VAF_", variable,"/",variable,"_discrete.png", sep = ""), dpi = 300, width = 15, height = 15, units = "in")
  
  
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
    
    # sub_mut <- subset(sub_mut, sub_mut$Sample != "PD8204a")
    
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
    
    ggsave(filename = paste("~/Desktop/MetaAML_results/Data/Figures/", "VAF_", variable,"/discrete/",variable, "_vs_", gene,"_VAF.png", sep = ""), dpi = 300, width = 5, height = 5, units = "in") 
    
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


# order the factors
var2_adj_list$Variable = factor(var2_adj_list$Variable, levels=c('WBC','Hemoglobin','Platelet','LDH', 'BM_blast_percent', 'PB_blast_percent', 'Age'))

# plot the discrete results
p = ggplot(var2_adj_list, aes(x = Mutated_Gene, y = effect_size, label = var2_adj_list$p_value)) +
  geom_hline(yintercept=0, linetype = "dashed", color = "black") +
  theme_cowplot() +
  # geom_text(aes(Mutated_Gene, var2_adj_list$CI_upper), hjust = 0, nudge_y = 0.1) +
  geom_pointrange(size = 0.75, stat = "identity", 
                  # shape = 21, fill = "white",
                  aes(x = Mutated_Gene, ymin = CI_lower, ymax = CI_upper, y = effect_size, color = sig)) +
  scale_color_manual(values = c("q < 0.2" = "#FDE725FF", "q < 0.1" = "#1F968BFF", "q < 0.05" = "#482677FF", "q > 0.2" = "grey")) +
  ylab("Effect Size (high VAF vs. low VAF)")+
  xlab("")+
  theme(legend.position = "right", legend.title = element_blank(),
        axis.title.y=element_blank()) +
  coord_flip() +
  xlim(rev(levels(factor(var2_adj_list$Mutated_Gene))))

p + facet_grid(. ~ Variable) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/VAF_clinical_features_discrete.png", dpi = 300, width = 15, height = 5, units = "in") 





# static vaf cutoff ####
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
                            risk.table = F,
                            risk.table.title = "",
                            risk.table.fontsize = 4,
                            risk.table.height = .3,
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
    png(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/",gene,"survival_by_median_VAF.png", sep = ""), res = 300, width = 3.5, height = 3.5, units = "in")
    
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
                            risk.table = F,
                            risk.table.title = "",
                            risk.table.fontsize = 4,
                            risk.table.height = .3,
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
    png(filename = paste("~/Desktop/MetaAML_results/Figure_3/Supplimental/",gene,"survival_by_30_VAF.png", sep = ""), res = 300, width = 3.5, height = 3.5, units = "in")
    
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
# temp_final_hr_order_15$categories = as.factor(temp_final_hr_order_15$categories)
# temp_final_hr_order_15$HR = as.numeric(temp_final_hr_order_15$HR)

temp_final = subset(temp_final, temp_final$Upper_CI != "Inf")

# temp_final$categories <- reorder(temp_final$categories, temp_final$HR)

temp_final$p_text = NA
temp_final$q_text = NA
for(i in 1:nrow(temp_final)){
  if(temp_final$log_rank_p[i] < 0.05){
    temp_final$p_text[i] = temp_final$log_rank_p[i]
    temp_final$q_text[i] = temp_final$fdr[i]
  }
}
# temp_final$q_text = NA
# for(i in 1:nrow(temp_final)){
#   if(temp_final$fdr[i] < 0.2){
#     temp_final$q_text[i] = temp_final$fdr[i]
#   }
# }

temp_final$p_text = as.numeric(temp_final$p_text)
temp_final$q_text = as.numeric(temp_final$q_text)

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
  scale_color_manual(values = c("0" = "darkgrey", "1" = "darkred"))+
  ylab("Hazard Ratio")+
  theme_cowplot() +
  ylim(0,14) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip() 

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/gene_vaf_discrete_hr_forest_plot_de_novo_30.pdf", dpi = 300, width = 5, height = 7.5, units = "in")

write.csv(temp_final, "~/Desktop/MetaAML_results/Data/Tables/static_vaf_threshold_survival.csv")




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
# temp_final_hr_order_15$categories = as.factor(temp_final_hr_order_15$categories)
# temp_final_hr_order_15$HR = as.numeric(temp_final_hr_order_15$HR)

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
    temp_final_median$p_text[i] = paste("p = ", paste(temp_final_median$p_text[i]))
  }
}
  
temp_final_median$p_q_text = paste("p = ", temp_final_median$p_text, "; q =", temp_final_median$q_text, sep = "")

for(i in 1:nrow(temp_final_median)){
  if(temp_final_median$log_rank_p[i] > 0.05){
    temp_final_median$p_q_text[i] = ""
  }
  if(temp_final_median$log_rank_p[i] > 0.05){
    temp_final_median$p_text[i] = ""
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
  scale_color_manual(values = c("0" = "darkgrey", "1" = "darkred"))+
  ylab("Hazard Ratio")+
  theme_cowplot() +
  ylim(0,18)+
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/gene_vaf_discrete_hr_forest_plot_de_novo_median.pdf", dpi = 300, width = 5, height = 7.5, units = "in")

write.csv(temp_final, "~/Desktop/MetaAML_results/Data/Tables/static_vaf_threshold_survival.csv")







# optimal vaf cutoff ####
dir.create("~/Desktop/MetaAML_results/Figure_3/optimal_vaf_thresholds")

# now define the optimal VAF cutoff for each gene
z <- 1
temp_list <- list()

results_list = list()
n=1

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
    
    png(filename = paste("~/Desktop/MetaAML_results/Figure_3/optimal_vaf_thresholds/", gene, "_vaf_logrank_plot.png"),res = 300, width = 5, height = 5, units = "in")
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
    
    if(n_distinct(final_data_matrix_2_sub2$vaf_threshold) > 1){
      # create the survival data object
      
      OS_fit <- survfit(Surv(Time_to_OS, Censor) ~ 1, data=final_data_matrix_2_sub2)
      OS_trt_fit <- survfit(Surv(Time_to_OS, Censor) ~ vaf_threshold, data=final_data_matrix_2_sub2, conf.type = "log-log")
      
      p <- surv_pvalue(OS_trt_fit)$pval
      
      # print(p)
      # store the results in a dataframe
      temp <- data.frame(matrix(NA, nrow = 1, ncol = 3))
      names(temp) <- c("Gene", "vaf_threshold", "p_value")
      
      temp[1,1] <- gene
      temp[1,2] <- threshold
      temp[1,3] <- p
      
      z <- z + 1
      temp_list[[z]] <- temp
      
      cohorts <- c("over", "under")
      
      title <- paste(threshold, "VAF\nthreshold",sep = " ")
      
      if(p <= 0.05){
        # plots the survival
        surv_plot <- ggsurvplot(OS_trt_fit,
                                data = final_data_matrix_2_sub2,
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
                                ggtheme = theme_cowplot())
                                # (plot.title = element_text(hjust = 0.5)))
        
        print(surv_plot)
        png(filename = paste("~/Desktop/MetaAML_results/Figure_3/optimal_vaf_thresholds/",gene,"survival_by_VAF.png"), res = 300, width = 5, height = 5, units = "in")
        
        surv_plot
        print(surv_plot)
        dev.off()
      }
      # summarize results from a Cox model
      final_data_matrix_2_sub2$vaf_threshold = ifelse(final_data_matrix_2_sub2$vaf_threshold == "over", 1,0)
      
      model <- coxph( Surv(Time_to_OS, Censor) ~ vaf_threshold,
                      data = final_data_matrix_2_sub2 )
      
      # extract the informative data from the survival model
      array_dat = summary(model)$conf.int[1:4]
      array_dat[5] = gene
      
      # extract the log-rank p-value for the individual comparisons
      array_dat[6] = summary(model)$sctest[3]
      array_dat = array_dat[-2]
      
      forest_plot_data <- data.frame("Gene" = array_dat[4], "HR" = array_dat[1], "Lower_CI" = array_dat[2], "Upper_CI" = array_dat[3], "log_rank_p" = array_dat[5])
      
      results_list[[n]] <- forest_plot_data
      n=n+1   
    }
  }
}

temp_final = as.data.frame(do.call(rbind, temp_list))
temp_final$q_value <- p.adjust(temp_final$q_value, method = "fdr")

write.csv(temp_final,  file = "~/Desktop/MetaAML_results/Figure_3/optimal_vaf_thresholds/optimal_vaf_threshold_for_survival_prediction.csv", row.names = F)


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
    temp_final_hr$p_text[i] = paste("p = ", paste(temp_final_hr$p_text[i]))
  }
}

temp_final_hr$p_q_text = paste(temp_final_hr$p_text, "; q =", temp_final_hr$q_text, sep = "")

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

ggplot(temp_final_hr, aes(x = reorder(Gene, -HR), y = HR, label = temp_final_hr$p_q_text)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_text(aes(Gene, Upper_CI), hjust = 0, nudge_y = 0.5) +
  geom_pointrange(size = 0.75, stat = "identity", shape = 19,
                  # fill = "white",
                  aes(x = Gene, ymin = Lower_CI, ymax = Upper_CI, y = HR, color = sig_color)) +
  scale_color_manual(values = c("0" = "#737373", "1" = "#1b7837", "2" = "#762a83"))+
  ylab("Hazard Ratio")+
  # ylim(0,11.5)+
  theme_cowplot() +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip()


ggsave(filename = "~/Desktop/MetaAML_results/Figure_3/optimal_vaf_thresholds/gene_vaf_optimal_vaf_thresholds_hr_forest_plot_de_novo.pdf", dpi = 300, width = 5, height = 7.5, units = "in")


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

# p values
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(final_data_matrix_2_sub)

library(RColorBrewer)


pdf(file = "~/Desktop/MetaAML_results/Figure_3/Supplimental/correlation_by_vaf.pdf", width = 7.5, height = 7.5)

corrplot(corr_mat, is.corr = F, type="upper", order="hclust",tl.col="black", outline = F, addgrid.col = "lightgrey",
         col = brewer.pal(n = 8, name = "RdBu"), p.mat = p.mat, diag=FALSE,insig = "label_sig",
         sig.level = c(.001, .01, .1), pch.cex = .9, pch.col = "black", na.label = "square", na.label.col = "white")
dev.off()


# cluster based on the VAFs ####
final_data_matrix_2_sub <- final_data_matrix_2

final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, mut_freq_gene > 50 & Subset == "de_novo")
final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, Gene != "MLL")

final_data_matrix_2_sub <- select(final_data_matrix_2_sub, Sample, Gene, VAF)
final_data_matrix_2_sub$VAF <- as.numeric(as.character(final_data_matrix_2_sub$VAF))
final_data_matrix_2_sub$VAF <- round(final_data_matrix_2_sub$VAF, 3)
# dup <- as.data.frame(duplicated(final_data_matrix_2_sub[c("Sample", "Gene")]))

# identify patients where the same gene is mutated twice and generate a duplicate patient to represent this additional mutation
for(i in 1:nrow(final_data_matrix_2_sub)){
  n <- as.character(final_data_matrix_2_sub$Sample[i])
  m <- as.character(final_data_matrix_2_sub$Gene[i])
  sub1 <- subset(final_data_matrix_2_sub, final_data_matrix_2_sub$Sample == n & final_data_matrix_2_sub$Gene == m)
  if(nrow(sub1) > 1){
    for(j in 2:nrow(sub1)){
      rn <- row_number(sub1)
      final_data_matrix_2_sub$Sample[i] <- paste(n, "_", rn, sep="")
    }
  }
}

final_data_matrix_2_sub <- dcast(final_data_matrix_2_sub, Sample ~ Gene, value.var="VAF")

rownames(final_data_matrix_2_sub) <- final_data_matrix_2_sub$Sample
final_data_matrix_2_sub$Sample <- NULL
final_data_matrix_2_sub[is.na(final_data_matrix_2_sub)] <- 0

final_data_matrix_2_sub = t(data.matrix(final_data_matrix_2_sub))

# cluster and viauslize the data
install.packages("pheatmap")
library("pheatmap")

p = pheatmap(final_data_matrix_2_sub, 
             border_color = NA,
             show_colnames = F,
             cluster_rows = T,
             cluster_cols = T
             # kmeans_k = 10
)

ggsave(p, filename = "~/Desktop/MetaAML_results/Figure_3/Supplimental/cluster_by_vaf.pdf", dpi = 300, width = 5, height = 5, units = "in")

