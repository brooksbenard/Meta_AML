# ========================================================================================================================================== #
# Figure_1.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 11/02/2021
# Description: This script will perform all analyses and generate plots for Figure 1 and related suppliments for Benard et al. "Clonal architecture predicts clinical outcomes and drug response in acute myeloid leukemia"
# ========================================================================================================================================== #
# Figure 1A ####
# UpSet plot of cohort ####
# MetaAML_data_upset_plot
# 
# This script plots the overlap of different data types for the meta analysis cohort
#
# load required packages
# Package names
packages <- c("ggplot2", "scales" , "readxl", "reshape2", "cowplot", "dplyr", "tidyr", "plyr", "UpSetR", "muhaz", "data.table", "ggpubr", "RCurl", "reshape", "survival", "survminer", "gsubfn")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

dir.create("~/Desktop/MetaAML_results/Figure_1")
setwd("~/Desktop/MetaAML_results/Figure_1")

load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_2 = final_data_matrix

# download annotation data from BeatAML website
# download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx")

# annotate RNA-seq and drug screened samples
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)
BeatAML_sample_data_types = dplyr::select(BeatAML_sample_data_types, LabId, PatientId, totalDrug, rnaSeq)
BeatAML_sample_data_types$rna_seq = ifelse(BeatAML_sample_data_types$rnaSeq == "y", 1,0)
BeatAML_sample_data_types$drug_screening = ifelse(BeatAML_sample_data_types$totalDrug == "y", 1,0)
BeatAML_sample_data_types = dplyr::select(BeatAML_sample_data_types, LabId, rna_seq, drug_screening)
colnames(BeatAML_sample_data_types)[1] = "Sample"

# TCGA expression data annotations
tcga_aml_rpkm=read.table("~/Desktop/MetaAML_results/raw_data/laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt", header = T)
tcga_aml_rpkm = as.data.frame(colnames(tcga_aml_rpkm))
tcga_aml_rpkm = as.data.frame(tcga_aml_rpkm[-1,])
tcga_aml_rpkm$rna_seq = 1
colnames(tcga_aml_rpkm)[1] = "Sample"
tcga_aml_rpkm$Sample = gsub("\\.", "-", tcga_aml_rpkm$Sample)

# bind BeatAML and TCGA annotation data together
bound = rbind.fill(tcga_aml_rpkm, BeatAML_sample_data_types)
bound[is.na(bound)] <- 0 

# format matrix of data types
# sets we want: de novo, secondary, relapse, VAF, survival
subset_data = dplyr::select(final_data_matrix_2, Sample, Gene, VAF, Subset, Cohort, Time_to_OS)

# add rna-seq annotation for TCGA and BeatAML datasets
subset_data = left_join(subset_data, bound, by = "Sample")
subset_data[is.na(subset_data)] <- 0 

# remove low frequency data types
subset_data = subset_data[subset_data$Subset %in% c('de_novo', 'secondary', 'relapse'),]

distinct_samples = as.data.frame(unique(subset_data$Sample))

upset_matrix <- data.frame(matrix(0, nrow = nrow(distinct_samples), ncol = 9))

names(upset_matrix) <- c("Sample", "De novo", "Secondary", "Relapse", "VAF", "Survival", "RNA", "Drug screening", "Cohort")

upset_matrix$Sample = distinct_samples$`unique(subset_data$Sample)`


for(i in 1:nrow(upset_matrix)){
  print(i)
  sub = subset(subset_data, Sample == upset_matrix$Sample[i])
  
  upset_matrix$Cohort[i] = as.character(sub[1,5])
  
  if("de_novo" %in% sub$Subset){
    upset_matrix$`De novo`[i] = 1
  }
  if("secondary" %in% sub$Subset){
    upset_matrix$Secondary[i] = 1
  }
  if("relapse" %in% sub$Subset){
    upset_matrix$Relapse[i] = 1
  }
  if(!is.na(sub$VAF)){
    upset_matrix$VAF[i] = 1
  }
  if(!is.na(sub$Time_to_OS[1])){
    upset_matrix$Survival[i] = 1
  }
  if(1 %in% sub$rna_seq){
    upset_matrix$RNA[i] = 1
  }
  if(1 %in% sub$drug_screening){
    upset_matrix$`Drug screening`[i] = 1
  }
}

upset_matrix$Sample = NULL
upset_matrix$Cohort = as.factor(upset_matrix$Cohort)

write.csv(upset_matrix, "~/Desktop/MetaAML_results/Figure_1/cohort_upset_plot.csv")

mtxcol <- data.frame(x=rep(1:14,each=7), 
                     color=rep(c("#252525","#252525","darkred", "#252525", "#252525", "#252525", "#252525"),each=14))


# plot the overlap of avaliable data types
p = upset(upset_matrix,
          nsets = 7,
          order.by = "freq",
          sets.bar.color = c("#252525","#252525","darkred", "#252525", "#252525", "#252525", "#252525"),
          shade.color = c("#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd"),
          shade.alpha = 1,
          queries = list(
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au", "Hirsch", "Huet", "Lindsley", "Parkin", "Wang", "Welch")), color = "#00A087FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au", "Hirsch", "Huet", "Lindsley", "Parkin", "Wang")), color = "#E64B35FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au", "Hirsch", "Huet", "Lindsley", "Parkin")), color = "#DC0000FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au", "Hirsch", "Huet", "Lindsley")), color = "#7AA6DCFF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au", "Hirsch", "Huet")), color = "#B09C85FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au", "Hirsch")), color = "#7E6148FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li", "Au")), color = "#4DBBD5FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif", "Li")), color = "#8491B4FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti", "Greif")), color = "#F39B7FFF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg", "Majeti")), color = "#868686FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA", "Garg")), color = "#3C5488FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner", "TCGA")), color = "#EFC000FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil", "Tyner")), color = "#0073C2FF", active = T),
            list(query = elements, 
                 params = list("Cohort", c("Papaemmanuil")), color = "#CD534CFF", active = T))
)
pdf(file = "~/Desktop/MetaAML_results/Figure_1/cohort_upset_plot.pdf", width = 5.5, height = 4)
p
dev.off()
# ========================================================================================================================================== #


# ========================================================================================================================================== #
# Figure 1B ####
# Oncoprint of genomic features ####
# MetaAML_oncoprint
packages <- c("ggplot2", "cowplot", "reshape2", "dplyr", "data.table", "ggsci", "magick", "devtools", "readr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# establish uniform colors for each cohort
cohort_colors = list("Tyner" = "#0073C2FF", 
                     "TCGA" = '#EFC000FF', 
                     "Majeti" = '#868686FF', 
                     "Papaemmanuil" = "#CD534CFF", 
                     "Lindsley" = '#7AA6DCFF', 
                     "Wang" = '#E64B35FF',
                     "Au" = '#4DBBD5FF',
                     "Welch" = '#00A087FF',
                     "Garg" = '#3C5488FF',
                     "Greif" = "#F39B7FFF",
                     "Li" = "#8491B4FF",
                     "Shlush" = "#91D1C2FF",
                     "Parkin" = "#DC0000FF", 
                     "Hirsch" = "#7E6148FF", 
                     "Huet" = "#B09C85FF")

subset_color = list("de_novo" = "#C16622FF", 
                    "secondary" = "#767676FF", 
                    "relapse" = "#800000FF", 
                    "other" = "#FFA319FF", 
                    "Remission" = "#8A9045FF",
                    "Residual" = "#155F83FF", 
                    "therapy" = "#8F3931FF", 
                    "MDS" = "#58593FFF")

final_data_matrix_3 <- final_data_matrix

for(i in 1:nrow(final_data_matrix_3)){
  if(final_data_matrix_3$variant_type[i] == "other"){
    final_data_matrix_3$variant_type[i] = "Unknown"
  }
}

# plot only genes mutated at least 100 times in the cohort
final_data_matrix_3 <- subset(final_data_matrix_3, final_data_matrix_3$mut_freq_gene >= 100)

c = as.numeric(n_distinct(final_data_matrix_3$Sample, final_data_matrix_3$Cohort, final_data_matrix_3$Age, final_data_matrix_3$Sex, final_data_matrix_3$Risk, final_data_matrix_3$Subset, final_data_matrix_3$Time_to_OS))
r = as.numeric(n_distinct(final_data_matrix_3$Gene))

temp_dat <- data.frame(matrix(NA, nrow = r, ncol = c))
names(temp_dat) <- unique(final_data_matrix_3$Sample)
rownames(temp_dat) <- unique(final_data_matrix_3$Gene)

final_data_matrix_3$variant_type <- as.character(final_data_matrix_3$variant_type)

for(i in 1:nrow(final_data_matrix_3)){
  
  print(i)
  
  pt <- final_data_matrix_3$Sample[i]
  gene <- final_data_matrix_3$Gene[i]
  var_type <- final_data_matrix_3$variant_type[i]
  
  k = as.numeric(which(colnames(temp_dat) == pt))
  l = as.numeric(which(rownames(temp_dat) == gene))
  
  temp_dat[l,k] <- var_type
}
temp_dat[is.na(temp_dat)] <- ""
temp_dat <- as.matrix(temp_dat)

# color for the variant type annotation
col = c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")

anno_df <- unique(data.frame(final_data_matrix_3$Sample, final_data_matrix_3$Cohort, final_data_matrix_3$Age, final_data_matrix_3$Sex, final_data_matrix_3$Risk, final_data_matrix_3$Subset, final_data_matrix_3$Time_to_OS))

colnames(anno_df) <- c("Sample", "Cohort", "Age", "Sex", "Risk", "Subset", "Survival")

for(i in 1:nrow(anno_df)){
  if(!is.na(anno_df$Survival[i])){
    anno_df$Survival[i] <- "Yes"
  }
}
anno_df$Survival[is.na(anno_df$Survival)] <- "No"

Sex = anno_df[, "Sex"]
Age = as.numeric(anno_df[, "Age"])
Cohort = anno_df[, "Cohort"]
Risk = anno_df[, "Risk"]
Subset = anno_df[, "Subset"]
Survival = anno_df[, "Survival"]

ha = HeatmapAnnotation(Sex = Sex, Cohort = Cohort, Risk = Risk, Subset = Subset, Survival = Survival,
                       col = list(Sex = c("Male" = "#6a51a3", 
                                          "Female" = "#43a2ca"),
                                  Survival = c("Yes" = "#252525",
                                               "No" = "#f0f0f0"),
                                  Risk = c("Adverse" = "#E64B35FF",
                                           "Intermediate" = "#8491B4FF",
                                           "Favorable" = "#00A087FF", 
                                           "Unknown" = "#767676FF"),
                                  Subset = c("de_novo" = "#C16622FF", 
                                             "secondary" = "#767676FF", 
                                             "relapse" = "#800000FF", 
                                             "other" = "#FFA319FF", 
                                             "Remission" = "#8A9045FF",
                                             "Residual" = "#155F83FF", 
                                             "therapy" = "#8F3931FF", 
                                             "MDS" = "#58593FFF"),
                                  Cohort = c("Tyner" = "#0073C2FF", 
                                             "TCGA" = '#EFC000FF', 
                                             "Majeti" = '#868686FF', 
                                             "Papaemmanuil" = "#CD534CFF", 
                                             "Lindsley" = '#7AA6DCFF', 
                                             "Wang" = '#E64B35FF',
                                             "Au" = '#4DBBD5FF',
                                             "Welch" = '#00A087FF',
                                             "Garg" = '#3C5488FF',
                                             "Greif" = "#F39B7FFF",
                                             "Li" = "#8491B4FF",
                                             "Shlush" = "#91D1C2FF",
                                             "Parkin" = "#DC0000FF", 
                                             "Hirsch" = "#7E6148FF", 
                                             "Huet" = "#B09C85FF")),
                       annotation_height = unit(c(3.5, 3.5, 3.5, 3.5, 3.5), "mm"), 
                       show_annotation_name = T, 
                       annotation_legend_param = list(Sex = list(title = "Sex"),
                                                      Survival = list(title = "Survival"),
                                                      Risk = list(title = "Risk"),
                                                      Subset = list(title = "Subset"),
                                                      Cohort = list(title = "Cohort"))
)


pdf()

fig1B = oncoPrint(temp_dat, 
                 col = col,
                 # top_annotation_height = 2,
                 row_names_side= "right",
                 # show_pct = T,
                 bottom_annotation = ha,
                 alter_fun_is_vectorized = FALSE,
                 alter_fun = list(
                   background = function(...) NULL,
                   
                   Deletion = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                             gp = gpar(fill = col["Deletion"], col = "#374E55FF")),
                   INDEL = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                          gp = gpar(fill = col["INDEL"], col = "#DF8F44FF")),
                   Insertion = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                              gp = gpar(fill = col["Insertion"], col = "#00A1D5FF")),
                   ITD = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                        gp = gpar(fill = col["ITD"], col = "#79AF97FF")),
                   SNV = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                        gp = gpar(fill = col["SNV"], col = "#B24745FF")),
                   Splicing = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                             gp = gpar(fill = col["Splicing"], col = "#6A6599FF")),
                   Unknown = function(x, y, w, h) grid.rect(x, y, w*0.5, h*0.9, 
                                                            gp = gpar(fill = col["Unknown"], col = "#80796BFF"))
                 ))

print(fig1B)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/cohort_oncoprint_100.pdf", res = 300, width = 15, height = 6.5, units = "in")

fig1B
print(fig1B)
dev.off()

# write out raw data for source data in manuscript
source_data = final_data_matrix %>%
  select(Sample, Gene, variant_type, Subset, specimenType, Cohort, Sex, Risk, Time_to_OS)

write_csv(source_data, "~/Desktop/MetaAML_results/Figure_1/cohort_oncoprint_100.csv")





#### Supplimental ####
# treatment histories ####
# find the treatment data for cohorts where avaliable. Then, homogenize all of the therapy labels into broad categories for induction therapy and potential bone marrow transplant
install.packages('alluvial'); library('alluvial')
# Papaemmanuil
Papaemmanuil_treatment = read.table("~/Desktop/MetaAML_results/raw_data/AML_knowledge_bank_data_clinical.txt", stringsAsFactors = F)
Papaemmanuil_treatment = Papaemmanuil_treatment %>%
  select(PDID, Study1, VPA, ATRA_arm, TPL_o, TPL_type, OS) 

Papaemmanuil_treatment$Induction = "ICE"

for(i in 1:nrow(Papaemmanuil_treatment)){
  if(Papaemmanuil_treatment$VPA[i] == 1 & Papaemmanuil_treatment$ATRA_arm[i] == 1){
    Papaemmanuil_treatment$Induction[i] = "ICE+VPA+ATRA"
  }
  if(Papaemmanuil_treatment$VPA[i] == 0 & Papaemmanuil_treatment$ATRA_arm[i] == 1){
    Papaemmanuil_treatment$Induction[i] = "ICE+ATRA"
  }
}

Papaemmanuil_treatment$TPL_type[is.na(Papaemmanuil_treatment$TPL_type)] <- "None"

names(Papaemmanuil_treatment)[c(6:7)] = c("Transplant_type", "Survival")

for(i in 1:nrow(Papaemmanuil_treatment)){
  if(Papaemmanuil_treatment$Transplant_type[i] %in% c("ALLO", "FREMD", "HAPLO", "TPL_(Spenderart_unbekannt)")){
    Papaemmanuil_treatment$Transplant_type[i] = "Allo"
  }
  if(Papaemmanuil_treatment$Transplant_type[i] == "AUTO"){
    Papaemmanuil_treatment$Transplant_type[i] = "Auto"
  } 
}

Papaemmanuil_treatment$Transplant = ifelse(Papaemmanuil_treatment$Transplant_type == "None", "No", "Yes")

Papaemmanuil_treatment$Survival = "Yes"
Papaemmanuil_treatment$Consolidation = ifelse(Papaemmanuil_treatment$Transplant == "Yes", "Transplant", "Unknown")

Papaemmanuil_treatment_simple = ddplyr::rename(count(Papaemmanuil_treatment, Induction, Transplant,Transplant_type, Consolidation, Survival), Freq = n)


Papaemmanuil_treatment_final <- Papaemmanuil_treatment %>% left_join(Papaemmanuil_treatment_simple, by=c("Induction", "Transplant", "Transplant_type", "Consolidation", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 4) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

Papaemmanuil_treatment_final %>%
  mutate(
    Cohort = "Papaemmanuil",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> Papaemmanuil_treatment_final

png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Papaemmanuil_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(Papaemmanuil_treatment_final[,1:5], freq=Papaemmanuil_treatment_final$Freq,
         col = Papaemmanuil_treatment_final$k,
         # border = TCGA_treatment_final$k,
         # hide = TCGA_treatment_final$Freq < 10 ,
         cex = 0.7
)
dev.off()


# Beat AML
Tyner_treatment = read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)
Tyner_treatment = Tyner_treatment %>%
  select(LabId, PatientId, cumulativeChemo, typeInductionTx, cumulativeTreatmentTypes, overallSurvival)

# simlify the survival coding
Tyner_treatment$Survival = "No"

for(i in 1:nrow(Tyner_treatment)){
  if(!is.na(Tyner_treatment$overallSurvival[i])){
    if(Tyner_treatment$overallSurvival[i] > 1){
      Tyner_treatment$Survival[i] = "Yes"
    }   
  }
}


# pull out the different types of induction therapy and simplify the type of transplant
Tyner_treatment$typeInductionTx[is.na(Tyner_treatment$typeInductionTx)] <- "Unknown"
Tyner_treatment$cumulativeTreatmentTypes[is.na(Tyner_treatment$cumulativeTreatmentTypes)] <- "Unknown"

# find the cases where the induction type appears to be standard chemotherapy
Tyner_treatment$Induction = "Other"
for(i in 1:nrow(Tyner_treatment)){
  if(Tyner_treatment$typeInductionTx[i] == "Unknown" & grep("Standard Chemotherapy", Tyner_treatment$cumulativeTreatmentTypes) & Tyner_treatment$cumulativeChemo[i] == "y"){
    Tyner_treatment$Induction[i] = "7+3"
  }  
  if(Tyner_treatment$typeInductionTx[i] == "Standard Chemotherapy"){
    Tyner_treatment$Induction[i] = "7+3"
  }
}

# find cases with a transplant
Tyner_treatment$Transplant = ifelse(grepl("Bone Marrow Transplant",Tyner_treatment$cumulativeTreatmentTypes), "Yes", "No")
# find cases with targeted therapy
Tyner_treatment$Targeted = ifelse(grepl("Targeted",Tyner_treatment$cumulativeTreatmentTypes), "Yes", "No")
# find cases with other therapy
Tyner_treatment$Other = ifelse(grepl("Other",Tyner_treatment$cumulativeTreatmentTypes), "Yes", "No")
# find cases with palliative/supportive care
Tyner_treatment$Palliative = ifelse(grepl("Supportive/Palliative Care",Tyner_treatment$cumulativeTreatmentTypes), "Yes", "No")

Tyner_treatment$Consolidation = "None"
for(i in 1:nrow(Tyner_treatment)){
  if(Tyner_treatment$Targeted[i] == "Yes" & Tyner_treatment$Palliative[i] == "Yes" & Tyner_treatment$Other[i] == "Yes"){
    Tyner_treatment$Consolidation[i] = "Targeted+other"
  }
  if(Tyner_treatment$Targeted[i] == "Yes" & Tyner_treatment$Palliative[i] == "Yes" & Tyner_treatment$Other[i] == "No"){
    Tyner_treatment$Consolidation[i] = "Targeted+other"
  }
  if(Tyner_treatment$Targeted[i] == "Yes" & Tyner_treatment$Palliative[i] == "No" & Tyner_treatment$Other[i] == "Yes"){
    Tyner_treatment$Consolidation[i] = "Targeted+other"
  }
  if(Tyner_treatment$Targeted[i] == "Yes" & Tyner_treatment$Palliative[i] == "No" & Tyner_treatment$Other[i] == "No"){
    Tyner_treatment$Consolidation[i] = "Targeted"
  }
  if(Tyner_treatment$Targeted[i] == "No" & Tyner_treatment$Palliative[i] == "Yes" & Tyner_treatment$Other[i] == "No"){
    Tyner_treatment$Consolidation[i] = "Palliative"
  }
  if(Tyner_treatment$Targeted[i] == "No" & Tyner_treatment$Palliative[i] == "No" & Tyner_treatment$Other[i] == "Yes" & Tyner_treatment$Transplant[i] != "Yes"){
    Tyner_treatment$Consolidation[i] = "Other"
  }
  if(Tyner_treatment$Transplant[i] == "Yes"){
    Tyner_treatment$Consolidation[i] = "Transplant"
  }
}


Tyner_treatment$Transplant_type = ifelse(Tyner_treatment$Transplant == "Yes", "Unknown", "None")


Tyner_treatment_simple = dplyr::rename(count(Tyner_treatment, Induction, Consolidation, Transplant,Transplant_type, Survival), Freq = n)

Tyner_treatment_final <- Tyner_treatment %>% left_join(Tyner_treatment_simple, by=c("Induction", "Transplant", "Consolidation", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

Tyner_treatment_final %>%
  mutate(
    Cohort = "Tyner",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> Tyner_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/BeatAML_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(Tyner_treatment_final[,1:5], freq=Tyner_treatment_final$Freq,
         col = Tyner_treatment_final$k,
         # border = TCGA_treatment_final$k,
         # hide = TCGA_treatment_final$Freq < 10 ,
         cex = 0.7
)
dev.off()




# TCGA
TCGA_treatment = read.table("~/Desktop/MetaAML_results/raw_data/data_clinical_patient.txt", header = T, stringsAsFactors = F, sep = "\t")
TCGA_treatment = TCGA_treatment %>%
  select(PATIENT_ID, INDUCTION, TRANSPLANT_TYPE, OS_MONTHS)

# simlify the survival coding
TCGA_treatment$Survival = "No"

for(i in 1:nrow(TCGA_treatment)){
  if(!is.na(TCGA_treatment$OS_MONTHS[i])){
    if(TCGA_treatment$OS_MONTHS[i] >= 0){
      TCGA_treatment$Survival[i] = "Yes"
    }   
  }
}

# test case 
for(i in 1:nrow(TCGA_treatment)){
  if(TCGA_treatment$TRANSPLANT_TYPE[i] %in% c("Auto", "Auto, MUD", "Auto, sib Allo")){
    TCGA_treatment$Transplant_type[i] = "Auto"
  }
  if(TCGA_treatment$TRANSPLANT_TYPE[i] %in% c("MUD", "MUD, Auto, MUD", "MUD, MUD", "MUD, MUD, MUD")){
    TCGA_treatment$Transplant_type[i] = "Allo"
  }
  if(TCGA_treatment$TRANSPLANT_TYPE[i] %in% c("sib Allo", "sib Allo, sib Allo")){
    TCGA_treatment$Transplant_type[i] = "Allo"
  }
  if(TCGA_treatment$TRANSPLANT_TYPE[i] == 0){
    TCGA_treatment$Transplant_type[i] = "None"
  }
}

for(i in 1:nrow(TCGA_treatment)){
  if(TCGA_treatment$INDUCTION[i] %in% c("7+3+ATRA", "+3+3+PSC", "7+3, IT", "7+3+Genasense", "7+3, dauna", "7+3+dauno", "7+3+study drug", "7+3+3, gleevec", "7+3+AMD", "7+3+3, then 5+2+2", "7+3+3+PSC", "7+4+ATRA", "Decitabine then 7+3", "Revlmd then Decitbne,7+3,5+2")){
    TCGA_treatment$Induction[i] = "7+3 +"
  }
  if(TCGA_treatment$INDUCTION[i] %in% c("Hydrea, ATRA started", "Hydrea & Idarubicin", "hydrea, didnt get addl chemo")){
    TCGA_treatment$Induction[i] = "Other"
  }
  if(TCGA_treatment$INDUCTION[i] %in% c("Decitabine", "LBH/Decitabine", "Revlimid", "CLAM", "Cytarabine only", "Azacitidine", "low dose Ara C", "no treatment")){
    TCGA_treatment$Induction[i] = "Other"
  }
  if(TCGA_treatment$INDUCTION[i] %in% c("7+3", "7+3+3")){
    TCGA_treatment$Induction[i] = "7+3"
  }
}

TCGA_treatment$Transplant = ifelse(TCGA_treatment$TRANSPLANT_TYPE == 0, "No", "Yes")
TCGA_treatment$Consolidation = ifelse(TCGA_treatment$Transplant == "Yes", "Transplant", "Unknown")

TCGA_treatment_simple = dplyr::rename(count(TCGA_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

TCGA_treatment_final <- TCGA_treatment %>% left_join(TCGA_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

TCGA_treatment_final %>%
  mutate(
    Cohort = "TCGA",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> TCGA_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/TCGA_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(TCGA_treatment_final[,1:5], freq=TCGA_treatment_final$Freq,
         col = TCGA_treatment_final$k,
         # ,
         cex = 0.7
)
dev.off()





# Stanford
Majeti_treatment = read_excel("~/Desktop/Majeti_Lab/Data/Old_stuff/Combined_mutation_occurence/151102_Clinical Outcomes of All Patients without_PHI copy.xlsx", sheet = 3)
Majeti_treatment = Majeti_treatment %>%
  select(`Patient SU Number`, T1, Induction, Consolidation, Transplant, `Duration from diagnosis to D1`) %>%
  subset(`Patient SU Number` %in% published_samples$Stanford_ID)

Majeti_treatment$Induction[is.na(Majeti_treatment$Induction)] <- "Unknown"
Majeti_treatment$Consolidation[is.na(Majeti_treatment$Consolidation)] <- "Unknown"

Majeti_treatment$Survival = "Yes"
Majeti_treatment$Transplant_type = "None"

for(i in 1:nrow(Majeti_treatment)){
  if(Majeti_treatment$Transplant[i] == 1){
    Majeti_treatment$Transplant_type[i] = "Unknown"
  }
  if(Majeti_treatment$Induction[i] == "GCLAC"){
    Majeti_treatment$Induction[i] = "Other"
    Majeti_treatment$Consolidation[i] = "Other"
  }
  if(Majeti_treatment$Consolidation[i] == "GCLAC"){
    Majeti_treatment$Consolidation[i] = "Other"
  }
  if(Majeti_treatment$Induction[i] %in% c("3+4", "AraC/Etoposide/Daunorubicin", "CPX-351")){
    Majeti_treatment$Induction[i] = "Other"
  }
  if(Majeti_treatment$Consolidation[i] %in% c("2+3", "2+5", "5+2", "3+5", "MEC")){
    Majeti_treatment$Consolidation[i] = "Other"
  }
}

Majeti_treatment_simple = dplyr::rename(count(Majeti_treatment, Induction, Consolidation, Transplant,Transplant_type, Survival), Freq = n)

Majeti_treatment_final <- Majeti_treatment %>% left_join(Majeti_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

Majeti_treatment_final %>%
  mutate(
    Cohort = "Majeti",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> Majeti_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Stanford_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(Majeti_treatment_final[,1:4], freq=Majeti_treatment_final$Freq,
         col = Majeti_treatment_final$k,
         
         cex = 0.7
)
dev.off()



# Lindsley
Lindsley_treatment = read_excel("~/Desktop/MetaAML_results/Data/Blood_2014_Lindsley_additional_clinical.xlsx")

Lindsley_treatment$Induction = "Other"
Lindsley_treatment$Consolidation = "Unknown"
Lindsley_treatment$Transplant = "Unknown"
Lindsley_treatment$Transplant_type = "Unknown"
Lindsley_treatment$Survival = "Yes"

Lindsley_treatment_simple = dplyr::rename(count(Lindsley_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

Lindsley_treatment_final <- Lindsley_treatment %>% left_join(Lindsley_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

Lindsley_treatment_final %>%
  mutate(
    Cohort = "Lindsley",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> Lindsley_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Lindsley_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(Lindsley_treatment_final[,1:4], freq=Lindsley_treatment_final$Freq,
         col = Lindsley_treatment_final$k,
         
         cex = 0.7
)
dev.off()



# BMC_2016_Au 
Au <- read_excel("~/Desktop/MetaAML_results/Data/BMC_2016_Au.xlsx") %>%
  select("No.")

names(Au) = c("Patient")
Au$Induction = "Unknown"
Au$Consolidation = "Unknown"
Au$Transplant = "Unknown"
Au$Transplant_type = "Unknown"
Au$Survival = "No"

Au_simple = dplyr::rename(count(Au, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

Au_final <- Au %>% left_join(Au_simple, by=c("Induction", "Transplant", "Consolidation", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()


# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

Au_final %>%
  mutate(
    Cohort = "Au",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> Au_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Au_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(Au_final[,1:5], freq=Au_final$Freq,
         col = Au_final$k,
         
         cex = 0.7
)
dev.off()


# Oncotarget_2016_Wang 
wang <- read_excel("~/Desktop/MetaAML_results/Data/Oncotarget_2016_Wang.xlsx")
wang <- wang[-1,]
wang <- wang %>%
  select(`Sample ID`) %>%
  unique()

wang$Induction = "Unknown"
wang$Consolidation = "Unknown"
wang$Transplant = "Unknown"
wang$Transplant_type = "Unknown"
wang$Survival = "No"

wang_simple = dplyr::rename(count(wang, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

wang_final <- wang %>% left_join(wang_simple, by=c("Induction", "Transplant", "Consolidation", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()


# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

wang_final %>%
  mutate(
    Cohort = "Wang",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> wang_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Au_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(wang_final[,1:5], freq=wang_final$Freq,
         col = wang_final$k,
         
         cex = 0.7
)
dev.off()


# Hirsch
hirsch_clinical = read_excel("~/Desktop/MetaAML_results/Data/Hirsch_2016_Nat_Comm_clinical_data.xlsx", sheet = 1)

hirsch_clinical = hirsch_clinical %>%
  select(UPN, `intensive chemotherapy`, `Allogeneic BMT in first CR`, `time from diagnosis (days)`) %>%
  unique() 

hirsch_clinical$Induction = ifelse(hirsch_clinical$`intensive chemotherapy` == "yes", "7+3", "Other")
hirsch_clinical$Transplant = ifelse(grepl("yes", hirsch_clinical$`Allogeneic BMT in first CR`) == T, "Yes", "No")
hirsch_clinical$Transplant_type = ifelse(grepl("Yes", hirsch_clinical$Transplant) == T, "Allo", "None")
hirsch_clinical$Survival = ifelse(hirsch_clinical$`time from diagnosis (days)` > 0, "Yes", "No")

hirsch_clinical$Consolidation = ifelse(hirsch_clinical$Transplant == "Yes", "Transplant", "Unknown")


hirsch_clinical_simple = dplyr::rename(count(hirsch_clinical, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

hirsch_clinical_final <- hirsch_clinical %>% left_join(hirsch_clinical_simple, by=c("Induction", "Transplant", "Consolidation", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()


# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

hirsch_clinical_final %>%
  mutate(
    Cohort = "Hirsch",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> hirsch_clinical_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Hirsch_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(hirsch_clinical_final[,1:5], freq=hirsch_clinical_final$Freq,
         col = hirsch_clinical_final$k,
         
         cex = 0.7
)
dev.off()




# Welch
welch_treatment = read_excel("~/Desktop/MetaAML_results/Data/NEJM_2016_Welch_clinical.xlsx")
welch_treatment = welch_treatment[-117,]
welch_treatment = welch_treatment %>%
  select(UPN, `Prior therapy`, `Transplant (Yes=1, No=0)`, `Survival (days)`)

names(welch_treatment) = c("UPN", "Induction", "Transplant", "Survival")

welch_treatment$Induction[is.na(welch_treatment$Induction)] <- "Unknown"
welch_treatment$Transplant = ifelse(welch_treatment$Transplant == 1, "Yes", "No")

for(i in 1:nrow(welch_treatment)){
  if(welch_treatment$Transplant[i] == "No"){
    welch_treatment$Transplant_type[i] = "None"
  }
  if(welch_treatment$Transplant[i] == "Yes"){
    welch_treatment$Transplant_type[i] = "Unknown"
  }
}
welch_treatment$Survival = "Yes"
welch_treatment$Consolidation = ifelse(welch_treatment$Transplant == "Yes", "Transplant", "Unknown")

'%ni%' <- Negate('%in%')

for(i in 1:nrow(welch_treatment)){
  if(welch_treatment$Induction[i] %in% c("7+3, HiDAC, desatinib maintenance", "7+3, HiDAC, oral Clofarabine x 5 cycles", "7+3. MitoEC", "7+3. HiDAC. CLAG", "7+3. HiDAC. FLAG. NK cell infusion", "7+3. HIDAC. NK cell infusion. Decitabine 5 days x 1 cycle")){
    welch_treatment$Induction[i] = "7+3"
  }
  if(welch_treatment$Induction[i] %ni% c("7+3, HiDAC, desatinib maintenance", "7+3, HiDAC, oral Clofarabine x 5 cycles", "7+3. MitoEC", "7+3. HiDAC. CLAG", "7+3. HiDAC. FLAG. NK cell infusion", "7+3. HIDAC. NK cell infusion. Decitabine 5 days x 1 cycle", "Unknown", "7+3")){
    welch_treatment$Induction[i] = "Other"
  }
}

welch_treatment_simple = dplyr::rename(count(welch_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

welch_treatment_final <- welch_treatment %>% left_join(welch_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

welch_treatment_final %>%
  mutate(
    Cohort = "Welch",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> welch_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Welch_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(welch_treatment_final[,1:4], freq=welch_treatment_final$Freq,
         col = welch_treatment_final$k,
         
         cex = 0.7
)
dev.off()


# Garg
garg_treatment_1 = read_excel("~/Desktop/MetaAML_results/Data/Garg_2015_blood_clinical_discovery_cohort.xlsx") %>%
  select(`Sample ID`, `Induction Chemotherapy`, `OS\r\n(Mon)`)

garg_treatment_2 = read_excel("~/Desktop/MetaAML_results/Data/Garg_2015_blood_clinical_targeted_cohort.xlsx") %>%
  select(`Sample ID`, `Induction Chemotherapy`, `OS\r\n(Mon)`)


garg_treatment = rbind(garg_treatment_1, garg_treatment_2)
garg_treatment$Transplant = "Unknown"
names(garg_treatment) = c("Patient", "Induction", "OS", "Transplant")

garg_treatment$Consolidation = "Unknown"
garg_treatment$Transplant_type = "Unknown"

garg_treatment$Survival = ifelse(garg_treatment$OS == "NA", "No", "Yes")

garg_treatment$Induction = "7+3"

garg_treatment_simple = dplyr::rename(count(garg_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

garg_treatment_final <- garg_treatment %>% left_join(garg_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()


# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

garg_treatment_final %>%
  mutate(
    Cohort = "Garg",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> garg_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Garg_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(garg_treatment_final[,1:5], freq=garg_treatment_final$Freq,
         col = garg_treatment_final$k,
         
         cex = 0.7
)
dev.off()



# Huet
huet_treatment = read_excel("~/Desktop/MetaAML_results/Data/Huet_2018_Blood.xlsx")
colnames(huet_treatment) = huet_treatment[1,]
huet_treatment = huet_treatment[-1,]
huet_treatment = huet_treatment %>%
  select(`Patient ID`, `Induction regimen`, `Time to final outcome (months)`)

names(huet_treatment) = c("Patient", "Induction", "OS")
huet_treatment$Consolidation = "Unknown"
huet_treatment$Transplant = "Unknown"
huet_treatment$Transplant_type = "Unknown"

for(i in 1:nrow(huet_treatment)){
  if(huet_treatment$Induction[i] == "3+7"){
    huet_treatment$Induction[i] = "7+3"
  }
  if(huet_treatment$Induction[i] == "3+7+ATRA"){
    huet_treatment$Induction[i] = "7+3 +"
  }
  if(huet_treatment$Induction[i] == "sequential"){
    huet_treatment$Induction[i] = "Other"
  }
}

huet_treatment$Survival = ifelse(huet_treatment$OS > 0, "Yes", "No")

huet_treatment_simple = dplyr::rename(count(huet_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

huet_treatment_final <- huet_treatment %>% left_join(huet_treatment_simple, by=c("Induction", "Transplant", "Consolidation", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()


# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

huet_treatment_final %>%
  mutate(
    Cohort = "Huet",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> huet_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Huet_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(huet_treatment_final[,1:5], freq=huet_treatment_final$Freq,
         col = huet_treatment_final$k,
         
         cex = 0.7
)
dev.off()





# Hirsch
hirsch_treatment = read_excel("~/Desktop/MetaAML_results/Data/Hirsch_2016_Nat_Comm_clinical_data.xlsx", sheet = 1) %>%
  select(UPN, `intensive chemotherapy`, `Allogeneic BMT in first CR`, `time from diagnosis (days)`)
names(hirsch_treatment) = c("Patient", "Induction_simple", "TPL", "OS")

hirsch_treatment$Induction = ifelse(hirsch_treatment$Induction_simple == "yes", "Intensive chemo", "None")
hirsch_treatment$Transplant = ifelse(grepl("yes", hirsch_treatment$TPL) == T, "Yes", "No")
hirsch_treatment$Transplant_type = ifelse(hirsch_treatment$Transplant == "Yes", "Allo", "None")

hirsch_treatment$Consolidation = ifelse(hirsch_treatment$Transplant == "Yes", "Transplent", "Unknown")

hirsch_treatment$Survival = ifelse(hirsch_treatment$OS > 0, "Yes", "No")

hirsch_treatment_simple = dplyr::rename(count(hirsch_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival), Freq = n)

hirsch_treatment_final <- hirsch_treatment %>% left_join(hirsch_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Freq) %>%
  subset(Freq > 2) %>%
  unique()

# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba")
pal <- c("#404040", "#2b83ba")

hirsch_treatment_final %>%
  mutate(
    Cohort = "Hirsch",
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> hirsch_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Hirsch_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(hirsch_treatment_final[,1:5], freq=hirsch_treatment_final$Freq,
         col = hirsch_treatment_final$k,
         
         cex = 0.7
)
dev.off()


# Azizi
"%ni%" <- Negate("%in%") 

multiple_treatment = read_excel("~/Desktop/MetaAML_results/Data/azizi_diagnosis_relapse_aml_clinical.xlsx") %>%
  select(Case, Cohort, Induction, Consolidation, Overall_Survival, Cohort) %>%
  subset(Cohort == "Grief")

multiple_treatment$Induction[is.na(multiple_treatment$Induction)] <- "Other"
multiple_treatment$Consolidation[is.na(multiple_treatment$Consolidation)] <- "Unknown"

multiple_treatment$Transplant = "No"
multiple_treatment$Transplant_type = "None"

for(i in 1:nrow(multiple_treatment)){
  if(grepl("Allo", multiple_treatment$Consolidation[i]) == T){
    multiple_treatment$Transplant[i] = "Yes"
    multiple_treatment$Transplant_type[i] = "Allo"
  }
  if(grepl("Auto", multiple_treatment$Consolidation[i]) == T){
    multiple_treatment$Transplant[i] = "Yes"
    multiple_treatment$Transplant_type[i] = "Auto"
  }
  if(multiple_treatment$Induction[i] %in% c("3+7","Cytarabine+Daunorubicin", "Daunarubicin+Cytarabine", "Daunarubicin+Cytarabine; Daunarubicin+Cytarabine")){
    multiple_treatment$Induction[i] = "7+3"
  }
  if(multiple_treatment$Induction[i] %in% c("3+7 + cyclosporine", "3+7, 2+5", "3+7 + gemtuzumab ozogamicin", "Cytarabine+Daunorubicin; Mitoxantrone+Cytarabine+Vincristine", "Cytarabine+Daunorubicin+Etoposide", "Cytarabine+Daunorubicin+Mitoxantrone", "Daunarubicin+Cytarabine (+/-bevacizumab); Cytarabine (+/-bevacizumab)", "Daunorubicin+Cytarabine+Etoposide; Daunorubicin+Cytarabine+ Etoposide")){
    multiple_treatment$Induction[i] = "7+3 +"
  }
  if(multiple_treatment$Induction[i] %in% c("Thioguanine+Cytarabine+Daunorubicin+Mitoxantrone", "Thioguanine+Cytarabine+Daunorubicin", "Mitozantrone+Cytarabine", "Idarubicin+Cytarabine", "Idarubicin+Cytarabine+Etoposide")){
    multiple_treatment$Induction[i] = "Other"
  }
  if(multiple_treatment$Consolidation[i] %in% c("Idarubicin+Cytarabine+Etoposide", "Amsacrine+Cytarabine; AutoSCT", "Amsacrine+Cytarabine; AlloSCT", "Amsacrine+Cytarabine", "Cytarabine")){
    multiple_treatment$Consolidation[i] = "Other"
  }
}


multiple_treatment$Overall_Survival[is.na(multiple_treatment$Overall_Survival)] <- "0"
multiple_treatment$Survival = ifelse(multiple_treatment$Overall_Survival > 0, "Yes", "No")

multiple_treatment_simple = dplyr::rename(count(multiple_treatment, Induction, Consolidation, Transplant, Transplant_type, Survival, Cohort), Freq = n)

multiple_treatment_final <- multiple_treatment %>% left_join(multiple_treatment_simple, by=c("Induction", "Consolidation", "Transplant", "Transplant_type", "Survival", "Cohort")) %>%
  select(Induction, Consolidation, Transplant, Transplant_type, Survival, Cohort, Freq) %>%
  subset(Freq > 2) %>%
  unique()


# pal <- c("#ca0020", "#f4a582","#bababa", "#404040", "#2b83ba", "darkred")
pal <- c("#404040", "#2b83ba")

multiple_treatment_final %>%
  mutate(
    ss = paste(Transplant),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> multiple_treatment_final


png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Azizi_treatment.png", res = 300, width = 10, height = 7.5, units = "in")
alluvial(multiple_treatment_final[,c(6,1:5)], freq=multiple_treatment_final$Freq,
         col = multiple_treatment_final$k,
         
         cex = 0.7
)
dev.off()



# all cohorts together
all_cohorts = dplyr::bind_rows(Papaemmanuil_treatment_final, Tyner_treatment_final, TCGA_treatment_final, Majeti_treatment_final,garg_treatment_final,huet_treatment_final,multiple_treatment_final,welch_treatment_final, Lindsley_treatment_final, Au_final, hirsch_clinical_final, wang_final)

for(i in 1:nrow(all_cohorts)){
  if(all_cohorts$Transplant_type[i] == "No"){
    all_cohorts$Transplant_type[i] = "None"
  }  
}

# "Welch" = '#00A087FF',
# "Wang" = '#E64B35FF',
# "Tyner" = "#0073C2FF", 
# "TCGA" = '#EFC000FF', 
# "Papaemmanuil" = "#CD534CFF", 
# "Majeti" = '#868686FF', 
# "Lindsley" = '#7AA6DCFF', 
# "Huet" = "#B09C85FF"
# "Greif" = "#F39B7FFF",
# "Garg" = '#3C5488FF',
# "Hirsch" = "#7E6148FF", 
# "Au" = '#4DBBD5FF',


pal <- c("#4DBBD5FF", "#7E6148FF", "#3C5488FF", "#F39B7FFF", "#B09C85FF",  "#7AA6DCFF", "#868686FF", "#CD534CFF", "#EFC000FF", "#0073C2FF", "#E64B35FF", "#00A087FF")

all_cohorts %>%
  mutate(
    ss = paste(Cohort),
    k = pal[ match(ss, sort(unique(ss))) ]
  ) -> all_cohorts_final

png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/Combined_treatment_all.png", res = 300, width = 30, height = 12, units = "in")
alluvial(all_cohorts_final[,c(7,1:5)], freq=all_cohorts_final$Freq,
         col = all_cohorts_final$k,
         # border = "darkgrey",
         cex = 1
)
dev.off()


# stacked barplot
Papaemmanuil_treatment_bar = Papaemmanuil_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation"))  %>%
  mutate(Cohort = "Papaemmanuil")
Tyner_treatment_bar = Tyner_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation"))  %>%
  mutate(Cohort = "Tyner")
TCGA_treatment_bar = TCGA_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation"))  %>%
  mutate(Cohort = "TCGA")
Majeti_treatment_bar = Majeti_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Majeti") 
garg_treatment_bar = garg_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Garg")
huet_treatment_bar = huet_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Huet")
multiple_treatment_bar = multiple_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Greif")
welch_treatment_bar = welch_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Welch")
Lindsley_treatment_bar = Lindsley_treatment %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Lindsley")
Au_bar = Au %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Au")
hirsch_treatment_bar =  hirsch_clinical %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Hirsch")
wang_bar =  wang %>%
  select(1, Transplant_type, Survival, Induction, Transplant, Consolidation) %>%
  magrittr::set_names(c("Sample", "Transplant_type", "Survival", "Induction", "Transplant", "Consolidation")) %>%
  mutate(Cohort = "Wang")

# all cohorts together
all_cohorts_bar = rbind(Papaemmanuil_treatment_bar, Tyner_treatment_bar, TCGA_treatment_bar, Majeti_treatment_bar,garg_treatment_bar,huet_treatment_bar,multiple_treatment_bar,welch_treatment_bar, Lindsley_treatment_bar, Au_bar, hirsch_treatment_bar, wang_bar)

all_cohorts_bar = all_cohorts_bar %>%
  subset(Cohort %in% c("Tyner" , 
                       "TCGA" , 
                       "Majeti", 
                       "Papaemmanuil", 
                       "Lindsley", 
                       "Wang" ,
                       "Au" ,
                       "Welch",
                       "Garg" ,
                       "Greif" ,
                       "Hirsch" , 
                       "Huet" )) 

for(i in 1:nrow(all_cohorts_bar)){
  if(all_cohorts_bar$Induction[i] %in% c("3+5", "Cytarabine+Mitoxantrone")){
    all_cohorts_bar$Induction[i] = "Other"
  }
}

# induction

a = ggplot(data=all_cohorts_bar, aes(x=forcats::fct_infreq(Induction), fill = Cohort)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values=c("Tyner" = "#0073C2FF", 
                             "TCGA" = '#EFC000FF', 
                             "Majeti" = '#868686FF', 
                             "Papaemmanuil" = "#CD534CFF", 
                             "Lindsley" = '#7AA6DCFF', 
                             "Wang" = '#E64B35FF',
                             "Au" = '#4DBBD5FF',
                             "Welch" = '#00A087FF',
                             "Garg" = '#3C5488FF',
                             "Greif" = "#F39B7FFF",
                             "Hirsch" = "#7E6148FF", 
                             "Huet" = "#B09C85FF")) +
  theme_cowplot() +
  ylab("# of patients")+
  xlab(NULL) +
  ggtitle("Induction") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/induction_distribution.pdf", dpi = 300, width = 7.5, height = 5, units = "in")


b = ggplot(data=all_cohorts_bar, aes(x=forcats::fct_infreq(Consolidation), fill = Cohort)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values=c("Tyner" = "#0073C2FF", 
                             "TCGA" = '#EFC000FF', 
                             "Majeti" = '#868686FF', 
                             "Papaemmanuil" = "#CD534CFF", 
                             "Lindsley" = '#7AA6DCFF', 
                             "Wang" = '#E64B35FF',
                             "Au" = '#4DBBD5FF',
                             "Welch" = '#00A087FF',
                             "Garg" = '#3C5488FF',
                             "Greif" = "#F39B7FFF",
                             "Hirsch" = "#7E6148FF", 
                             "Huet" = "#B09C85FF")) +
  theme_cowplot() +
  ylab("# of patients")+
  xlab(NULL) +
  ggtitle("Consolidation") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/consolidation_distribution.pdf", dpi = 300, width = 7.5, height = 5, units = "in")


c = ggplot(data=all_cohorts_bar, aes(x=forcats::fct_infreq(Transplant), fill = Cohort)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values=c("Tyner" = "#0073C2FF", 
                             "TCGA" = '#EFC000FF', 
                             "Majeti" = '#868686FF', 
                             "Papaemmanuil" = "#CD534CFF", 
                             "Lindsley" = '#7AA6DCFF', 
                             "Wang" = '#E64B35FF',
                             "Au" = '#4DBBD5FF',
                             "Welch" = '#00A087FF',
                             "Garg" = '#3C5488FF',
                             "Greif" = "#F39B7FFF",
                             "Hirsch" = "#7E6148FF", 
                             "Huet" = "#B09C85FF")) +
  theme_cowplot() +
  ylab("# of patients")+
  xlab(NULL) +
  ggtitle("Transplant") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5))

ggsave(filename = "~/DesktopMetaAML_results/Figure_1/Supplimental/transplant_distribution.pdf", dpi = 300, width = 7.5, height = 5, units = "in")

d = ggplot(data=all_cohorts_bar, aes(x=forcats::fct_infreq(Transplant_type), fill = Cohort)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values=c("Tyner" = "#0073C2FF", 
                             "TCGA" = '#EFC000FF', 
                             "Majeti" = '#868686FF', 
                             "Papaemmanuil" = "#CD534CFF", 
                             "Lindsley" = '#7AA6DCFF', 
                             "Wang" = '#E64B35FF',
                             "Au" = '#4DBBD5FF',
                             "Welch" = '#00A087FF',
                             "Garg" = '#3C5488FF',
                             "Greif" = "#F39B7FFF",
                             "Hirsch" = "#7E6148FF", 
                             "Huet" = "#B09C85FF")) +
  theme_cowplot() +
  ylab("# of patients")+
  xlab(NULL) +
  ggtitle("Transplant Type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(
    strip.background = element_rect(colour="black", fill="white", 
                                    size=1.5, linetype="solid"),
    plot.title = element_text(hjust = 0.5))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/transplant_type_distribution.pdf", dpi = 300, width = 7.5, height = 5, units = "in")

ggarrange(a, b, c, d, 
          ncol = 4,
          legend = "right",
          align = "hv",
          # align = "v",
          common.legend = TRUE)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/treatment_distributions.pdf", dpi = 300, width = 20, height = 5, units = "in")



if (!require('ggridges')) install.packages('ggridges'); library('ggridges')

dir.create("~/Desktop/MetaAML_results/Figure_1/Supplimental")

# data
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# establish uniform colors for each cohort
cohort_colors = c("Tyner" = "#0073C2FF", 
                  "TCGA" = '#EFC000FF', 
                  "Majeti" = '#868686FF', 
                  "Papaemmanuil" = "#CD534CFF", 
                  "Lindsley" = '#7AA6DCFF', 
                  "Wang" = '#E64B35FF',
                  "Au" = '#4DBBD5FF',
                  "Welch" = '#00A087FF',
                  "Garg" = '#3C5488FF',
                  "Greif" = "#F39B7FFF",
                  "Li" = "#8491B4FF",
                  "Shlush" = "#91D1C2FF",
                  "Parkin" = "#DC0000FF", 
                  "Hirsch" = "#7E6148FF", 
                  "Huet" = "#B09C85FF")

# Age ####
sub=na.omit(as.data.frame(distinct(final_data_matrix, Sample, Age, Cohort)))
sub$Age=as.numeric(sub$Age)

a1=ggplot(sub, aes(x = Age, y = as.factor(Cohort), fill = Cohort)) +
  theme_cowplot(font_size = 10) +
  geom_density_ridges(
    jittered_points = F, position = "raincloud",
    alpha = 1, scale = 1,
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_y_discrete(expand = c(0,0)) +
  ylab(label = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values = cohort_colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/age_by_cohort.png", dpi = 300, width = 3, height = 4, units = "in")

# Gender ####
sub=na.omit(as.data.frame(distinct(final_data_matrix, Sample, Sex, Cohort)))

sub = sub %>% group_by(Cohort, Sex) %>% tally()

b1=ggplot(sub, aes(fill=Sex, y=n, x=Cohort)) + 
  geom_bar(position = "fill", stat="identity",) +
  scale_fill_manual(values = c("Male" = "#6a51a3", "Female" = "#43a2ca")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  ylab(label = NULL) +
  xlab(label = NULL)  +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/sex_by_cohort.png", dpi = 300, width = 5, height = 3, units = "in")

# Risk ####
sub=na.omit(as.data.frame(distinct(final_data_matrix, Sample, Risk, Cohort)))

sub = sub %>% group_by(Cohort, Risk) %>% tally()

c1=ggplot(sub, aes(fill=Risk, y=n, x=Cohort)) + 
  geom_bar(position = "fill", stat="identity",) +
  scale_fill_manual(values = c("Adverse" = "#E64B35FF",
                               "Intermediate" = "#8491B4FF",
                               "Favorable" = "#00A087FF", 
                               "Unknown" = "#767676FF")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   # size = 12, 
                                   hjust = 1)) +
  ylab(label = NULL) +
  xlab(label = NULL)  +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/risk_by_cohort.png", dpi = 300, width = 5, height = 3, units = "in")


# add age, sex, and risk to the same grid
ggarrange(a1,              
          ggarrange(b1, c1, nrow = 2), 
          nrow = 1
)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/age_sex_risk.pdf", dpi = 300, width = 6, height = 3.25, units = "in")




# Survival ####
# de novo
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
sub = na.omit(as.data.frame(distinct(final_data_matrix, Sample, Cohort, Time_to_OS, Censor, Subset)))

sub = subset(sub, Subset == "de_novo")

sub$Time_to_OS <- (sub$Time_to_OS/365)
sub$OS <- with(sub, Surv(Time_to_OS, Censor == 1))
OS <- survfit(OS ~ Cohort, data = sub, conf.type = "log-log")

surv_plot <- ggsurvplot(OS,
                        data = sub,
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
                        risk.table.height = .3,
                        risk.table.y.text = T,
                        break.time.by = 5,
                        risk.table.pos = c("out"),
                        # palette = cohort_colors,
                        xlab = "Years",
                        ylim = c(0, 1.0),
                        ylab =  "Survival Probability",
                        font.main = c(15, "plain", "#252525"),
                        pval.size = 4,
                        font.x = c(12, "plain", "#252525"),
                        font.y =  c(12, "plain", "#252525"),
                        font.legend = c(12, "plain"),
                        font.tickslab = c(12, "plain", "#252525"),
                        legend.labs = c("Garg", "Hirsch", "Huet", "Majeti", "Papaemmanuil", "TCGA", "Tyner", "Welch"),
                        legend.title = "Cohort",
                        legend = "right",
                        ggtheme = theme_cowplot(font_size = 10))

print(surv_plot)
png(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/de_novo_survival.png", res = 300, width = 5, height = 3, units = "in")

surv_plot
print(surv_plot)
dev.off()


# forrest plot ####
library("survival")
library("survminer")
# individual mutation's HRs
sub = final_data_matrix %>%
  select(Sample, Gene, variant_type, Time_to_OS, Censor, Subset, mut_freq_gene) %>%
  subset(Subset == "de_novo" & mut_freq_gene > 75)
sub$Time_to_OS <- (sub$Time_to_OS/365)

# make sure that the FLT3 symbols are annotated well
for(i in 1:nrow(sub)){
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "ITD"){
    sub$Gene[i] <- "FLT3_ITD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "SNV"){
    sub$Gene[i] <- "FLT3_TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "Deletion"){
    sub$Gene[i] <- "FLT3_TKD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "INDEL"){
    sub$Gene[i] <- "FLT3_ITD"
  }
  if(sub$Gene[i] == "FLT3" & sub$variant_type[i] == "other"){
    sub$Gene[i] <- "FLT3_TKD"
  }
}

colnames(sub) = c("Sample", "Gene", "variant_type", "time", "status", "Subset", "mut_freq_gene")

sub = sub %>%
  select(Sample, Gene, time, status) %>%
  unique() %>%
  na.omit() 

mut_hr = data.frame(matrix(NA, nrow = length(unique(sub$Sample)), ncol = 30))
colnames(mut_hr) = c("Sample", "time", "status","NPM1","DNMT3A","SRSF2","TET2","KRAS","SF3B1","RUNX1","IDH1","NRAS","PTPN11","KIT","STAG2","IDH2","WT1","TP53","U2AF1","EZH2", "GATA2","BCOR","PHF6","CBL","ASXL1","CEBPA","RAD21", "MLL", "FLT3_ITD",  "FLT3_TKD")

mut_hr[,c(1:3)] = unique(sub[,c(1,3:4)])

for(i in 4:ncol(mut_hr)){
  gene = colnames(mut_hr[i])
  for(j in 1:nrow(mut_hr)){
    sample_gene = subset(sub, Sample == mut_hr$Sample[j])
    if(gene %in% sample_gene$Gene){
      mut_hr[j,i] = 1
    } else {
      mut_hr[j,i] = 0
    }
  }
}

mut_hr$status = as.numeric(as.character(mut_hr$status))
mut_hr = as.data.frame(mut_hr)

# set covariates for survival analysis
covariates <- c("NPM1","DNMT3A","SRSF2","TET2","KRAS","SF3B1","RUNX1","IDH1","NRAS","PTPN11","KIT","STAG2","IDH2","WT1","TP53","U2AF1","EZH2", "GATA2","BCOR","PHF6","CBL","ASXL1","CEBPA","RAD21", "MLL", "FLT3_ITD",  "FLT3_TKD")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste0('Surv(time, status)~', x, sep = "")))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mut_hr)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         # HR <- paste0(HR, " (", 
                         #              HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR", "lower_95", "upper_95", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res = as.data.frame(res)

res$gene = rownames(res)


res$gene <- factor(res$gene, levels = res$gene[order(res$HR)])
# res$log_rank_p = round(res$p.value, 3)
res$p_text = NA

for(i in 1:nrow(res)){
  if(res$p.value[i] <= 0.05){
    res$p_text[i] = paste("p = ", res$p.value[i], sep = "")
  }
  # if(res$log_rank_p[i] == 0){
  #   res$p_text[i] = paste("p < 0.001")
  # }
  if(res$log_rank_p[i] > 0.05){
    res$p_text[i] = ""
  }
}

# add functional category to the mutations for visualization purposes
res$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3_ITD", "FLT3_TKD", "JAK2")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33", "MLL")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(res)){
  if(res$gene[i] %in% DNA_methylation){
    res$mutation_category[i] <- "DNA Methylation"
  }
  if(res$gene[i] %in% Chromatin_cohesin){
    res$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(res$gene[i] %in% RTK_RAS_Signaling){
    res$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(res$gene[i] %in% Splicing){
    res$mutation_category[i] <- "Splicing"
  }
  if(res$gene[i] %in% Transcription){
    res$mutation_category[i] <- "Transcription"
  }
  if(res$gene[i] == "NPM1"){
    res$mutation_category[i] <- "NPM1"
  }
  if(res$gene[i] %in% Tumor_suppressors){
    res$mutation_category[i] <- "Tumor suppressors"
  }
}

# extract the p-value and hazard ratio for the individual interactions
# p = round(as.numeric(summary(model)$sctest[3]), 5)
# 
# # p = ifelse(p < 0.001, paste0("p < 0.001"), paste("p =", p))
# 
# hr = paste("HR = ", round(as.numeric(forest$HR), 2), " (", round(as.numeric(forest$lower_95), 2), "-", round(as.numeric(forest$upper_95), 2), ")", sep = "")
# 
# p_hr = paste(p, "; ", hr, sep = "")

# color coded by mutation category plot
k1 = ggplot(res, aes(x = reorder(gene, -HR), y = HR, label = res$p_text)) +
  geom_hline(yintercept=.5, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  geom_hline(yintercept=2, linetype="dashed", color = "#d9d9d9") +
  geom_hline(yintercept=3, linetype="dashed", color = "#d9d9d9") +
  geom_text(aes(gene, upper_95), hjust = 0, nudge_y = 0.35, size = 3) +
  theme_cowplot(font_size = 10) +
  ylim(0,5.25) +
  geom_pointrange(size = .75, stat = "identity", shape = 15, 
                  aes(x = gene, ymin = lower_95, ymax = upper_95, y = HR, color = mutation_category)) +
  scale_color_manual(name = "", values = c("DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF", "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF", "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF")) +
  ylab("Hazard Ratio")+
  xlab(NULL) +
  theme(legend.position = c(0.6,0.25),
        axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/individual_gene_HR_forrest_plot_de_novo.pdf", dpi = 300, width = 9, height = 5, units = "in")


# VAF per cohort ####
sub = na.omit(as.data.frame(distinct(final_data_matrix, Sample, Gene, VAF, Cohort, Subset)))

l1 = ggplot(sub, aes(x = VAF, y = as.factor(Cohort), fill = Cohort)) +
  theme_cowplot(font_size = 10) +
  geom_density_ridges(
    jittered_points = F, 
    position = "raincloud",
    alpha = 1, scale = 1,
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_y_discrete(expand = c(0,0)) +
  ylab(label = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values = cohort_colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/vaf_by_cohort.png", dpi = 300, width = 3, height = 4, units = "in")



# VAF per gene ####
sub <- subset(final_data_matrix, mut_freq_gene >= 75)
sub = na.omit(as.data.frame(distinct(sub, Sample, Gene, VAF, Cohort, Subset)))

m1 = ggplot(sub, aes(x = VAF, y = reorder(Gene, desc(Gene)))) +
  geom_density_ridges(
    jittered_points = F, 
    position = "raincloud",
    alpha = 1, scale = 1,
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_y_discrete(expand = c(0,0)) +
  xlim(0,100) +
  ylab(label = NULL) +
  theme(legend.position = "none") +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/vaf_by_gene.png", dpi = 300, width = 3, height = 4, units = "in")


# panel of features
# add age, sex, and risk to the same grid
ggarrange(k1,l1,m1,
          nrow = 1,
          ncol = 3,
          align = "hv",
          widths = c(1,.5,.5)
)

ggsave(filename = "~/Desktop/MetaAML_results/Figure_1/Supplimental/gene_hr_vaf_panel.pdf", dpi = 300, width = 15, height = 4.5, units = "in")
