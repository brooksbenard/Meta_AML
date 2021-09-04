# ========================================================================================================================================== #
# Figure_5.R
# Author : Brooks Benard and Logan Leak, bbenard@stanford.edu; lleak@stanford.edu
# Date: 08/23/2021
# Description: This script performes clonal analysis using PyClone and ClonEvol on aggregate TCGA and Beat AML samples
# This script will download, process, and generate the restults as seen in Figure 5 (and associated suppliments) of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #


# SECTIONS
# 1. Packages, libraries, and working directory
# 2. Reading in BeatAML dataset and keeping only desired columns
# 3. Implementing exclusion criteria on BeatAML
# 4. Creating columns for PyClone Analysis for BeatAML
# 5. Reading in TCGA dataset and keeping only desired columns
# 6. Implementing exclusion criteria on TCGA data
# 7. Creating columns for PyClone Analysis for TCGA data
# 8. Merging BeatAML and TCGA datasets
# 9. Exporting each sample into a unique data table in a unique folder in the pyclone_results_3 folder
# 10. Printing the script to be pasted into the terminal for PyClone
# 11. Run pasted script in terminal to get PyClone results for each patient
# 12. Reading in the data from the PyClone results and merging with original data frame
# 13. adding age, survival, and risk information of BeatAML to the PyClone data frame
# 14. adding age, survival, and risk information of TCGA to the PyClone data frame
# 15. Final merge of all data
# 16. Calculation of Shannon diversity scores
# 17. Removal of non-functional mutations
# 18. Graphs
# 19. Fisher's Exact Test
# 20. Preparing Figures
# 21. Quartiles


######################################################
# 1. Packages, libraries, and working directory   ####
######################################################
# Package names
packages <- c("ggplot2", "tidyverse", "viridis", "viridisLite", "RColorBrewer", "tydyr", "plyr", "dplyr", "cometExactTest", "discover", "stringr", "maditr", "reshape2", "data.table", "epitools", "corrplot", "plyr", "muhaz", "reshape", "survival", "survivalAnalysis", "survMisc", "survminer", "ggsci", "vegan", "ggrepel", "ggforce", "rstatix", "effsize", "psych", "maxstat", "RCurl", "ggpubr", "UpSetR", "cowplot", "readxl", "scales", "rlist", "ggplotify", "ggridges", "gridExtra", "magrittr", "stats", "gage", "fgsea", "GSA", "MAGeCKFlute", "devtools", "gridBase", "igraph", "packcircles")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Put last to avoid conflicts
install_github('hdng/clonevol')
library(clonevol)
if (!require('vegan')) install.packages('vegan'); library('vegan')


# Make directories
dir.create("~/Desktop/MetaAML_results/Figure_5")
dir.create("~/Desktop/MetaAML_results/Figure_5/PyClone")
setwd("~/Desktop/MetaAML_results/Figure_5/PyClone/")

#####################################################################
# 2. Reading in BeatAML dataset and keeping only desired columns ####
#####################################################################

## Download BeatAML dataset
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx")

## Reading in BeatAML mutation data file
BeatAML <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7, col_names = TRUE)

## Reading in BeatAML clinical data file
BeatAML_clinical <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5, col_names = TRUE)

## Merging datasets and consolidating columns

# Making labId column title in BeatAML uppercase so that it matches the clinical dataset
colnames(BeatAML)[which(colnames(BeatAML) == "labId")] <- "LabId"

# Full-joining the two datasets based on the labId
BeatAML <- merge(BeatAML, BeatAML_clinical, by = "LabId")

## Select only the columns that I need
BeatAML <- BeatAML %>%
  select(LabId, PatientId, consensus_sex, ageAtDiagnosis, chrom, pos_start, ref, alt, isDenovo, specimenType, total_reads, allele_reads, protein_position, amino_acids, all_consequences, ELN2017, vitalStatus, overallSurvival, isRelapse, symbol, genotyper, Karyotype)

######################################################
# 3. Implementing exclusion criteria on BeatAML   ####
######################################################

## Filter out patients less than 18 years old at diagnosis, patients that are not de novo, and patients that are relapsed

# exclude this for now and filter after (Brooks)

# BeatAML <- BeatAML %>%
#   filter(ageAtDiagnosis > 18) %>%
#   filter(isDenovo) %>%
#   filter(!isRelapse)

## Keeping only the hits that are found by both mutect and varscan. Also keeping all pindel hits

# Calculating VAF by dividing allele reads by total reads
BeatAML$VAF <- NA
BeatAML$VAF <- BeatAML$allele_reads/BeatAML$total_reads

# arrange by VAF in descending order so that we keep only the highest VAF from overlapping hits
BeatAML <- BeatAML %>%
  arrange(desc(VAF))

# create new data frame with the columns needed to determine mutation overlaps 
overlaps <- BeatAML %>%
  select(LabId, symbol, amino_acids, protein_position)

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
BeatAML <- BeatAML %>%
  arrange(LabId)

## keep only the labIds desired for PyClone. If there is more than one LabId for a patient, only keep the bone marrow LabId if present, otherwise keep all of the LabIds
unique_patient_ids <- unique(BeatAML$PatientId)
BeatAML$use <- NA
for (i in 1:length(unique_patient_ids)) {
  if (length(unique(BeatAML$LabId[which(BeatAML$PatientId == unique_patient_ids[i])])) > 1 && ("Bone Marrow Aspirate" %in% BeatAML$specimenType[which(BeatAML$PatientId == unique_patient_ids[i])])) {
    BeatAML$use[intersect(which(BeatAML$PatientId == unique_patient_ids[i]), which(BeatAML$specimenType == "Bone Marrow Aspirate"))] <- TRUE
    BeatAML$use[intersect(which(BeatAML$PatientId == unique_patient_ids[i]), which(BeatAML$specimenType != "Bone Marrow Aspirate"))] <- FALSE
    # if there is more than one LabId for a sample and Bone Marrow Aspirate is present, only use bone marrow aspirate data and put FALSE for the rest
  } else {
    BeatAML$use[which(BeatAML$PatientId == unique_patient_ids[i])] <- TRUE
    # else keep all of the LabIds for that patient
  }
}

# remove unwanted values
BeatAML <- BeatAML[which(BeatAML$use),]


###########################################################
# 4. Creating columns for PyClone Analysis for BeatAML ####
###########################################################

## Columns include:
# 1. mutation_id
# 2. ref_counts
# 3. var_counts
# 4. normal_cn
# 5. minor_cn
# 6. major_cn

## 1. mutation_id
# mutation id will be a concatenation of chrom, pos_start, ref, and alt
BeatAML$mutation_id <- NA
for (i in 1:nrow(BeatAML)) {
  BeatAML$mutation_id[i] <- str_c(BeatAML$chrom[i], BeatAML$pos_start[i], BeatAML$ref[i], BeatAML$alt[i], sep = "_")
}

## 2. ref_counts
# subtract allele_reads column from total_reads column
BeatAML$ref_counts <- NA
for (i in 1:nrow(BeatAML)) {
  BeatAML$ref_counts[i] <- BeatAML$total_reads[i] - BeatAML$allele_reads[i]
}

## 3. var_counts
# make a copy of allele reads called var_counts
BeatAML$var_counts <- NA
BeatAML$var_counts <- BeatAML$allele_reads

## 4. normal_cn
# if consensus_sex is male and chrom is X or Y, then it is 1. Otherwise it is 2
BeatAML$normal_cn <- NA
for (i in 1:nrow(BeatAML)) {
  if (!is.na(BeatAML$consensus_sex[i]) && !is.na(BeatAML$chrom[i])) {
    if (BeatAML$consensus_sex[i] == "Male" && (BeatAML$chrom[i] == "X" || BeatAML$chrom[i] == "Y")) {
      BeatAML$normal_cn[i] <- 1
    } else {
      BeatAML$normal_cn[i] <- 2
    }
  } else {
    BeatAML$normal_cn[i] <- 2
  }
}

#####################################
# perform copy number imputation using avaliable karyotype information

# download a file correlating the chromosome loci with gene ID
# http://uswest.ensembl.org/biomart/martview/a3043f537a26692111e9a7a94003ff68
all_genes <- read.table("~/Downloads/mart_export.txt", sep = "\t", header = T, fill = T,  stringsAsFactors = FALSE, quote = "")

all_genes$Karyotype_arm = gsub("\\..*","",all_genes$Karyotype.band)
all_genes$partial_annotation = paste("(",all_genes$Chromosome.scaffold.name,")","(",all_genes$Karyotype_arm,")", sep = "")
all_genes$full_annotation = paste("(",all_genes$Chromosome.scaffold.name,")","(",all_genes$Karyotype.band,")", sep = "")

# append to mutation file
names(all_genes)[1] = "symbol"

# append the chromosomal arm/band column to the mutation table
cohort_aggrigate = left_join(BeatAML, all_genes, by = "symbol")

## 5. minor_cn
# # Make equal to zero if no CNA information
cohort_aggrigate$minor_cn <- 0
# 
# ## 6. major_cn
# # Make equal to two if CNA information
cohort_aggrigate$major_cn <- 2

# create columns for new VAF based on copy number correction
cohort_aggrigate$VAF_CN_corrected = NA

# loop through all mutations and assign copy num ber status
for(i in 1:nrow(cohort_aggrigate)){
  
  # find whole chromosome gains or losses
  ch_gain = paste("\\+",cohort_aggrigate$Chromosome.scaffold.name[i], ",", sep = "")
  ch_loss = paste("\\-",cohort_aggrigate$Chromosome.scaffold.name[i], ",", sep = "")
  
  # correct for VAF based on copy number gains for that chromosome  
  if(grepl(ch_gain, cohort_aggrigate$Karyotype[i]) & cohort_aggrigate$VAF[i] > .35) {
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*0.75
    cohort_aggrigate$minor_cn[i] = 1
    cohort_aggrigate$major_cn[i] = 2
  }
  if(grepl(ch_gain, cohort_aggrigate$Karyotype[i]) & cohort_aggrigate$VAF[i] <= .35) {
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*1.5
    cohort_aggrigate$minor_cn[i] = 1
    cohort_aggrigate$major_cn[i] = 2
  }
  
  # correct for VAF based on copy number loss for that chromosome
  if(grepl(ch_loss, cohort_aggrigate$Karyotype[i])) {
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2
    cohort_aggrigate$minor_cn[i] = 0
    cohort_aggrigate$major_cn[i] = 1
  }
  
  # find focal chromosome gains or losses
  locus_gain = paste("add",cohort_aggrigate$partial_annotation[i],"|","add",cohort_aggrigate$full_annotation[i], sep = "")
  locus_loss = paste("del",cohort_aggrigate$partial_annotation[i],"|","del",cohort_aggrigate$full_annotation[i], sep = "")
  
  # correct for VAF based on broad copy number gains at the gene locus  
  if(grepl(locus_gain, cohort_aggrigate$Karyotype[i], fixed = T) & cohort_aggrigate$VAF[i] > .35) {
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*0.75
    cohort_aggrigate$minor_cn[i] = 1
    cohort_aggrigate$major_cn[i] = 2
    }
  if(grepl(locus_gain, cohort_aggrigate$Karyotype[i], fixed = T) & cohort_aggrigate$VAF[i] <= .35){
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]*1.5
    cohort_aggrigate$minor_cn[i] = 1
    cohort_aggrigate$major_cn[i] = 2
    }
  # correct for VAF based on broad copy number loss at the gene locus
  if(grepl(locus_loss, cohort_aggrigate$Karyotype[i], fixed = T)) {
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2
    cohort_aggrigate$minor_cn[i] = 0
    cohort_aggrigate$major_cn[i] = 1
    }
}

# TP53 is a unique case. Code manually for focal deletions
tp53_loss1 = "del(17)(p13)"
tp53_loss2 = "del(17)(p11.2p13)"

for(i in 1:nrow(cohort_aggrigate)){
  # if no copy number differenes detected, populate the raw VAF
  if(grepl(tp53_loss1, cohort_aggrigate$Karyotype[i]) & cohort_aggrigate$symbol[i] == "TP53") { 
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2
    cohort_aggrigate$minor_cn[i] = 0
    cohort_aggrigate$major_cn[i] = 1
    }
  
  if(grepl(tp53_loss2, cohort_aggrigate$Karyotype[i], fixed = T) & cohort_aggrigate$symbol[i] == "TP53") { 
    cohort_aggrigate$VAF_CN_corrected[i] = cohort_aggrigate$VAF[i]/2
    cohort_aggrigate$minor_cn[i] = 0
    cohort_aggrigate$major_cn[i] = 1
    }
}

BeatAML = cohort_aggrigate

######################################


##################################################################
# 5. Reading in TCGA dataset and keeping only desired columns ####
##################################################################

## Downloading data file
download.file("http://download.cbioportal.org/laml_tcga_pub.tar.gz", destfile = "~/Desktop/MetaAML_results/raw_data/laml_tcga_pub.tar.gz")

## unzip folder and extract mutation, clinical, and copy number alteration files
untar("~/Desktop/MetaAML_results/raw_data/laml_tcga_pub.tar.gz", files = c("data_mutations_extended.txt", "data_clinical_patient.txt", "data_CNA.txt"), exdir = "raw_data/")

## read in raw mutation file
TCGA <- arrange(read.table(file = '~/Desktop/MetaAML_results/raw_data/data_mutations_extended.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, quote = ""), Matched_Norm_Sample_Barcode)

## read in clinical data
TCGA_clinical <- arrange(read.table("~/Desktop/MetaAML_results/raw_data/data_clinical_patient.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE), PATIENT_ID)

## rename Matched_Norm_Sample_Barcode on TCGA dataframe so that they can be merged
TCGA$PATIENT_ID <- NA
TCGA$PATIENT_ID <- TCGA$Matched_Norm_Sample_Barcode

## merge datasets
TCGA <- merge(TCGA, TCGA_clinical, by = "PATIENT_ID")

## select only the columns I need
TCGA <- TCGA %>%
  select(PATIENT_ID, SEX, AGE, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, t_ref_count, t_alt_count, Variant_Classification, RISK_CYTO, OS_STATUS, OS_MONTHS, Hugo_Symbol)

## Change long names in Variant_Classification column
Missense <- c(858, 847, 409, 134, 184, 630, 1382, 35, 325, 1308, 1061, 829, 1353, 1352)
Silent <- c(351, 1378)

for (num in Missense) {
  TCGA$Variant_Classification[num] <- "Missense_Mutation"
}

for (num in Silent) {
  TCGA$Variant_Classification[num] <- "Silent"
}

######################################################
# 6. Implementing exclusion criteria on TCGA data ####
######################################################

## Filter out patients less than 18 years old at diagnosis
TCGA <- TCGA %>%
  filter(AGE > 18)

## filter out patients missing t_ref_count or t_alt_count
TCGA <- TCGA[union(which(!is.na(TCGA$t_alt_count)), which(!is.na(TCGA$t_ref_count))),]

#############################################################
# 7. Creating columns for PyClone Analysis for TCGA data ####
#############################################################

## Columns include:
# 1. mutation_id
# 2. ref_counts
# 3. var_counts
# 4. normal_cn
# 5. minor_cn
# 6. major_cn

## 1. mutation_id
# mutation_id will be a concatenation of Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2
TCGA$mutation_id <- NA

for (i in 1:nrow(TCGA)) {
  TCGA$mutation_id[i] <- str_c(TCGA$Chromosome[i], TCGA$Start_Position[i], TCGA$Reference_Allele[i], TCGA$Tumor_Seq_Allele2[i], sep = "_")
}

## 2. ref_counts
# rename column t_ref_count to ref_counts
colnames(TCGA)[which(colnames(TCGA) == "t_ref_count")] <- "ref_counts"

## 3. var_counts
# rename column t_alt_count to var_counts
colnames(TCGA)[which(colnames(TCGA) == "t_alt_count")] <- "var_counts"

## 4. normal_cn
# if SEX is male and Chromosome is X or Y, then normal_cn is 1. Otherwise, it is 2
TCGA$normal_cn <- NA
for (i in 1:nrow(TCGA)) {
  if (TCGA$SEX[i] == "Male" && (TCGA$Chromosome[i] == "X" || TCGA$Chromosome[i] == "Y")) {
    TCGA$normal_cn[i] <- 1
  } else {
    TCGA$normal_cn[i] <- 2
  }
}

## 5. minor_cn
TCGA$minor_cn <- 0

## 6. major_cn
# read in copy number alteration data frame
CNA <- read.table("~/Desktop/MetaAML_results/raw_data/data_CNA.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# must convert PATIENT_IDs from format TCGA-AB-2802 to TCGA.AB.2802 to match CNA data frame
TCGA$CNA_ID <- NA

for (i in 1:nrow(TCGA)) {
  current_id <- as.character(TCGA$PATIENT_ID[i])
  TCGA$CNA_ID[i] <- str_c(substr(current_id, 1, 4), substr(current_id, 6, 7), substr(current_id, 9, 12), "03", sep = ".")
}

# incorporate data from CNA data frame into major_cn calculation
TCGA$major_cn <- NA

for (i in 1:nrow(TCGA)) {
  current_id <- TCGA$CNA_ID[i]
  current_symbol <- TCGA$Hugo_Symbol[i]
  delta <- as.integer(CNA[which(CNA$Hugo_Symbol == current_symbol), which(colnames(CNA) == current_id)])
  if (length(delta) == 0) {
    delta <- 0
  }
  TCGA$major_cn[i] <- 2 + delta
}

###########################################
# 8. Merging BeatAML and TCGA datasets ####
###########################################

## rename BeatAML PatientId to PATIENT_ID so that all column names match
colnames(BeatAML)[which(colnames(BeatAML) == "LabId")] <- "PATIENT_ID"
colnames(TCGA)[which(colnames(TCGA) == "Hugo_Symbol")] <- "symbol"

# TCGA <- rename(TCGA, symbol = Hugo_Symbol)

# Stratify FLT3 into FLT3-ITD and FLT3-TKD
BeatAML <- BeatAML %>%
  mutate(symbol = ifelse(symbol == "FLT3" & genotyper == "pindel", "FLT3-ITD", symbol)) %>%
  mutate(symbol = ifelse(symbol == "FLT3" & genotyper != "pindel", "FLT3-TKD", symbol))

TCGA <- TCGA %>%
  mutate(symbol = ifelse(symbol == "FLT3" & PATIENT_ID == "TCGA-AB-2970", "FLT3-ITD", symbol)) %>%
  mutate(symbol = ifelse(symbol == "FLT3" & PATIENT_ID != "TCGA-AB-2970", "FLT3-TKD", symbol))

PyClone <- unique(rbind(select(BeatAML, PATIENT_ID, mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, symbol), select(TCGA, PATIENT_ID, mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, symbol)))


###################################################
# 9. Exporting each sample into a unique data     #
# table in a unique folder in the                 #
# pyclone_results_3 folder                        #
###################################################
dir.create("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders")

## creating a list of all of the future file names
all_files <- paste0(unique(PyClone$PATIENT_ID), ".tsv")

## creating a folder for each file within the PyClone_results_3 folder
for (i in 1:length(all_files)) {
  file <- all_files[i]
  if (substring(file, 1, 1) == "T") {
    dir.create(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", substring(file, 1, 12)))
  } else {
    dir.create(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", substring(file, 1, 8)))
  }
}

## creating a unique .tsv file for each patient and putting it in the appropriate folder
all_samples <- unique(PyClone$PATIENT_ID)

for (id in all_samples) {
  print(id)
  write.table(data.frame("mutation_id" = PyClone$mutation_id[which(PyClone$PATIENT_ID == id)], "ref_counts" = PyClone$ref_counts[which(PyClone$PATIENT_ID == id)], "var_counts" = PyClone$var_counts[which(PyClone$PATIENT_ID == id)], "normal_cn" = PyClone$normal_cn[which(PyClone$PATIENT_ID == id)], "minor_cn" = PyClone$minor_cn[which(PyClone$PATIENT_ID == id)], "major_cn" = PyClone$major_cn[which(PyClone$PATIENT_ID == id)]), paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", paste0(id), "/", paste0(id), ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)
}

###################################################
# 10. Printing the script to be pasted into the   #
# terminal for PyClone                            #
###################################################

result <- ""
command <- "pyclone run_analysis_pipeline --in_files "
for (i in 1:length(all_files)) {
  file <- all_files[i]
  if (substring(file, 1, 1) == "T") {
    result <- paste0(result, "cd ~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", substring(file, 1, 12), " ; ", command, all_files[i], " --working_dir ~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", substring(file, 1, 12), " ; ")
  } else {
    result <- paste0(result, "cd ~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", substring(file, 1, 8), " ; ", command, all_files[i], " --working_dir ~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", substring(file, 1, 8), " ; ")
  }
}

write(result, "~/Desktop/MetaAML_results/Figure_5/PyClone/pyclone_run_script.txt")


#############################################################################################
# 11. Run pasted script from 10 in your terminal to get PyClone results for each patient ####
#############################################################################################
# You can install PyClone using bioconda.
# 
# conda install pyclone -c bioconda -c conda-forge
# 
# This will install PyClone into your current conda environment. In some cases it may be better to create a separate conda environment for PyClone which be activated when needed. This avoids issues due to conflicting libraries. To create the environment execute the following command.
# 
# conda create -n pyclone -c bioconda -c conda-forge pyclone
# 
# Once the environment is created it can be activated using the following command.
# 
# conda activate pyclone
# 
# You can check that PyClone was installed correctly by running the following command which will show the help.
# 
# PyClone --help

############################################################################################
# 12. Reading in the data from the PyClone results and merging with original data frame ####
############################################################################################
## There will be three data frames from PyClone:
# 1. PyClone_clonality - has the number of clones and total number of mutations per patient
# 2. PyClone_granular - has the mutation_id, sample_id, cluster_id, cellular_prevalence,	cellular_prevalence_std,	and variant_allele_frequency (as determine by PyClone) for each sample
# 3. PyClone_cluster - has the sample_id, cluster_id, size, mean, and std for each clone

#####################
# PyClone_clonality #
#####################

# creating a data frame with the number of clones and total number of mutations per patient
PyClone_clonality <- data.frame("PATIENT_ID" = all_samples, stringsAsFactors = FALSE)
PyClone_clonality$number_of_clones <- NA
PyClone_clonality$total_mutations <- NA

# reading in data
for (id in all_samples) {
  print(id)
  if(id != "TCGA-AB-2820" & id !="TCGA-AB-2886" & id !="TCGA-AB-2937" & id !="TCGA-AB-2945" & id !="TCGA-AB-2966" & id !="TCGA-AB-3009"){
    setwd(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", paste0(id), "/tables"))
    if (file.exists("cluster.tsv")) {
      table <- read.table(file = "cluster.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      PyClone_clonality[which(PyClone_clonality$PATIENT_ID == id), which(colnames(PyClone_clonality) == "number_of_clones")] <- length(table$cluster_id)
      PyClone_clonality[which(PyClone_clonality$PATIENT_ID == id), which(colnames(PyClone_clonality) == "total_mutations")] <- sum(table$size)
    } else {
      PyClone_clonality[which(PyClone_clonality$PATIENT_ID == id), which(colnames(PyClone_clonality) == "number_of_clones")] <- 1
      PyClone_clonality[which(PyClone_clonality$PATIENT_ID == id), which(colnames(PyClone_clonality) == "total_mutations")] <- 1 
  }
  }
}

# set working directory back to PyClone_results_3
setwd("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders")

## merge the clonality data into PyClone data frame
PyClone <- full_join(PyClone, PyClone_clonality, by = "PATIENT_ID")


####################
# PyClone_granular #
####################
PyClone_granular <- data.frame(mutation_id = character(), sample_id = character(), cluster_id = integer(), cellular_prevalence = double(),	cellular_prevalence_std = double(),	variant_allele_frequency = double())

for (id in all_samples) {
  if(id != "TCGA-AB-2820" & id !="TCGA-AB-2886" & id !="TCGA-AB-2937" & id !="TCGA-AB-2945" & id !="TCGA-AB-2966" & id !="TCGA-AB-3009"){
  setwd(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", paste0(id), "/tables/"))
  }
  if (file.exists("loci.tsv")) {
    table <- read.table(file = "loci.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    PyClone_granular <- rbind(PyClone_granular, table)
  }
}

## merge PyClone_granular with PyClone
# rename sample_id column to PATIENT_ID for merge
colnames(PyClone_granular)[which(colnames(PyClone_granular) == "sample_id")] <- "PATIENT_ID"

PyClone <- full_join(PyClone, PyClone_granular, by = c("mutation_id", "PATIENT_ID"))

###################
# PyClone_cluster #
###################
PyClone_cluster <- data.frame(sample_id = character(), cluster_id = integer(), size = integer(), mean = double(), std = double())
for (id in all_samples) {
  if(id != "TCGA-AB-2820" & id !="TCGA-AB-2886" & id !="TCGA-AB-2937" & id !="TCGA-AB-2945" & id !="TCGA-AB-2966" & id !="TCGA-AB-3009"){
  setwd(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/PyClone_folders/", paste0(id), "/tables/"))
  }
  if (file.exists("cluster.tsv")) {
    table <- read.table(file = "cluster.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    PyClone_cluster <- rbind(PyClone_cluster, table)
  }
}
# renaming PyClone_cluster sample_id column name to PATIENT_ID
colnames(PyClone_cluster)[which(colnames(PyClone_cluster) == "sample_id")] <- "PATIENT_ID"

# merge with PyClone data
PyClone <- full_join(PyClone, PyClone_cluster, by = c("PATIENT_ID", "cluster_id"))


####################################################################################################
# 13. Preparing to add age, survival, and risk information of BeatAML to the PyClone data frame ####
####################################################################################################

# rename BeatAML age column name to AGE
colnames(BeatAML)[which(colnames(BeatAML) == "ageAtDiagnosis")] <- "AGE"

# rename BeatAML vitalStatus to Censor
colnames(BeatAML)[which(colnames(BeatAML) == "vitalStatus")] <- "Censor"

# overallSurvival does not need to be renamed

# rename BeatAML ELN2017 to Risk
colnames(BeatAML)[which(colnames(BeatAML) == "ELN2017")] <- "Risk"

#################################################################################################
# 14. Preparing to add age, survival, and risk information of TCGA to the PyClone data frame ####
#################################################################################################

# AGE column is correct

# rename TCGA OS_STATUS column to Censor
colnames(TCGA)[which(colnames(TCGA) == "OS_STATUS")] <- "Censor"
colnames(BeatAML)[which(colnames(BeatAML) == "vitalStatus")] <- "Censor"

# convert TCGA OS_MONTHS to days. Average number of days per month is 30.4375
TCGA$OS_MONTHS <- TCGA$OS_MONTHS*30.4375
# change TCGA OS_MONTHS column name to overallSurvival
colnames(TCGA)[which(colnames(TCGA) == "OS_MONTHS")] <- "overallSurvival"

# rename TCGA RISK_CYTO column to Risk
colnames(TCGA)[which(colnames(TCGA) == "RISK_CYTO")] <- "Risk"

##################################
# 15. Final merge of all data ####
##################################
Age_survival <- rbind(distinct(select(TCGA, PATIENT_ID, AGE, Censor, overallSurvival, Risk)), distinct(select(BeatAML, PATIENT_ID, AGE, Censor, overallSurvival, Risk)))
PyClone <- full_join(PyClone, Age_survival, by = "PATIENT_ID")

# make censor values consistent between BeatAML and TCGA (0 = alive, 1 = deceased)
for (i in 1:length(PyClone$PATIENT_ID)) {
  if (str_detect(PyClone$Censor[i], "(Alive)|(LIVING)")) {
    PyClone$Censor[i] <- 0
  } else if (str_detect(PyClone$Censor[i], "(Dead)|(DECEASED)")) {
    PyClone$Censor[i] <- 1
  } else {
    PyClone$Censor[i] <- NA
  }
}

write_tsv(PyClone, "~/Desktop/MetaAML_results/Figure_5/PyClone_summary_table.tsv")

PyClone = unique(PyClone)

##################################################
# 16. Calculation of Shannon diversity scores ####
##################################################
## Create new data frame with PATIENT_ID, cluster_id, number_of_clones, and mean
Shannon <- PyClone %>%
  select(PATIENT_ID, cluster_id, number_of_clones, mean) %>%
  distinct() %>%
  filter(number_of_clones > 1) %>%
  filter(!is.na(number_of_clones)) %>%
  filter(!is.na(mean))

## Calculate normalized mean for each cluster
Shannon$normalized_mean <- NA
for (i in 1:nrow(Shannon)) {
  Shannon$normalized_mean[i] <- Shannon$mean[i]/sum(Shannon$mean[which(Shannon$PATIENT_ID == Shannon$PATIENT_ID[i])])
}

## Calculate Shannon Diversity Index for each sample
Shannon$SDI <- NA
for (i in 1:nrow(Shannon)) {
  Shannon$SDI[i] <- diversity(Shannon$normalized_mean[which(Shannon$PATIENT_ID == Shannon$PATIENT_ID[i])], index = "shannon", MARGIN = 1, base = exp(1))
}


## merge into PyClone
PyClone <- Shannon %>%
  select(PATIENT_ID, SDI) %>%
  distinct() %>%
  full_join(PyClone, by = "PATIENT_ID")


symbols <- PyClone %>%
  select(symbol) %>%
  na.omit() %>%
  group_by(symbol) %>%
  summarise(frequency = length(symbol)) %>%
  arrange(desc(frequency))
print(symbols$symbol)


# PyClone = read_tsv("~/Desktop/Figure_5/PyClone/raw_data/summary_table.tsv")

##############################################
# 17. Removal of non-functional mutations ####
##############################################
## rename BeatAML all_consequences column name to Variant_Classification to match TCGA
colnames(BeatAML)[which(colnames(BeatAML) == "all_consequences")] <- "Variant_Classification"

## merge together BeatAML and TCGA mutation type data and merge into PyClone data frame
PyClone <- BeatAML %>%
  select(PATIENT_ID, mutation_id, Variant_Classification) %>%
  rbind(select(TCGA, PATIENT_ID, mutation_id, Variant_Classification)) %>%
  full_join(PyClone, by = c("PATIENT_ID", "mutation_id"))

# add annotations for cohort, subset, etc.
beat_aml_labels = BeatAML_clinical %>%
  select(LabId, PatientId, isRelapse, isDenovo, isTransformed)

beat_aml_labels$Cohort = "Tyner"

names(beat_aml_labels)[1] = "PATIENT_ID"

PyClone = left_join(PyClone, beat_aml_labels, by = "PATIENT_ID")

PyClone$Cohort = PyClone$Cohort %>% replace_na("TCGA")
PyClone$isDenovo = PyClone$isDenovo %>% replace_na("TRUE")
PyClone$isRelapse = PyClone$isRelapse %>% replace_na("FALSE")

for(i in 1:nrow(PyClone)){
  PyClone$PatientId[i] = PyClone$PatientId[i] %>% replace_na(PyClone$PATIENT_ID[i])  
}

PyClone$isDenovo = as.logical(PyClone$isDenovo)
PyClone$isRelapse = as.logical(PyClone$isRelapse)

PyClone_final = PyClone

## create vector of types of mutations to remove
# remove_these <- c("splice_donor_variant", "splice_acceptor_variant", "splice_donor_variant&non_coding_transcript_variant", "splice_acceptor_variant&intron_variant", "IGR", "Intron", "RNA", "Silent")

remove_these <- c("splice_donor_variant&non_coding_transcript_variant", "splice_acceptor_variant&intron_variant", "IGR", "Intron", "RNA", "Silent")

## remove any matching mutations
PyClone$Keep <- NA
for (i in 1:nrow(PyClone)) {
  if (PyClone$Variant_Classification[i] %in% remove_these) {
    PyClone$Keep[i] <- FALSE
  } else if (is.na(PyClone$Variant_Classification[i])) {
    PyClone$Keep[i] <- FALSE
  } else {
    PyClone$Keep[i] <- TRUE
  }
}
PyClone_with_all_mutations <- PyClone

## Remove false values
PyClone <- PyClone[which(PyClone$Keep),]

## Reassign total mutation number for each patient to only count the remaining functional mutations and also add a driver mutation count
for (i in 1:nrow(PyClone)) {
  PyClone$total_mutations[i] <- length(which(PyClone$PATIENT_ID == PyClone$PATIENT_ID[i]))
  PyClone$mutation_frequency[i] <- length(which(PyClone$symbol == PyClone$symbol[i]))
}

# calculate the number of driver mutations per patient
for (i in 1:nrow(PyClone)) {
  PyClone$n_driver_mutations[i] <- length(which(PyClone$PATIENT_ID == PyClone$PATIENT_ID[i] & PyClone$mutation_frequency >= 3))
}


symbols <- PyClone %>%
  select(symbol) %>%
  na.omit() %>%
  group_by(symbol) %>%
  summarise(frequency = length(symbol)) %>%
  arrange(desc(frequency))
print(symbols$symbol)



PyClone_final = PyClone
write_tsv(PyClone_final, "~/Desktop/MetaAML_results/Figure_5/PyClone_final.tsv")

PyClone = read_tsv("~/Desktop/MetaAML_results/Figure_5/PyClone_final.tsv")


# filter final results to only adult de novo patients and recurrent mutations
PyClone = PyClone %>%
  filter(AGE > 18) %>%
  filter(isDenovo) %>%
  filter(!isRelapse) %>%
  filter(n_driver_mutations > 0)

# for later clonal evolution analysis
PyClone_with_all_mutations = PyClone_with_all_mutations %>%
  filter(AGE > 18) %>%
  filter(isDenovo) %>%
  filter(!isRelapse)

PyClone$overallSurvival = PyClone$overallSurvival/365
PyClone_with_all_mutations$overallSurvival = PyClone_with_all_mutations$overallSurvival/365
# n_distinct(PyClone$PatientId)
# number of patients = 329

#################
# 18. Graphs ####
################# 
#### PyClone Clustering ####
# Figure 5B
# PyClone analysis and viz

# load packages
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('pheatmap')) install.packages('pheatmap'); library('pheatmap')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

# load data
df = read_tsv("~/Desktop/MetaAML_results/Figure_5/PyClone_final.tsv")

symbols <- df %>%
  select(symbol) %>%
  na.omit() %>%
  group_by(symbol) %>%
  summarise(frequency = length(symbol)) %>%
  arrange(desc(frequency)) %>%
  filter(frequency > 10) %>%
  as.data.frame()

# filter to only the most frequencly mutated genes
# genes_list = c("SPEN",  "CROCC",  "NRAS", "SMC3", "WT1",  "CBL",  "ETV6", "KRAS", "PTPN11", "FLT3", "FLT3-ITD", "FLT3-TKD", "IDH2", "TP53", "NF1",  "SRSF2",  "SETBP1", "CEBPA", "DNMT3A", "SF3B1",  "IDH1",   "ASXL1" , "RUNX1" , "U2AF1" , "EP300" , "GATA2" , "KIT" ,   "TET2"  , "NPM1",   "EZH2" ,  "RAD21" , "MYC" ,   "JAK2",   "ZRSR2" , "BCOR"  , "SMC1A",  "STAG2" , "BCORL1" ,"PHF6" ,  "NOTCH1" ,"MLL2" ,  "MLL")

df = subset(df, df$symbol %in% symbols$symbol)
df = subset(df, df$cellular_prevalence != "NA")

# create unique identifiers for the clusters per patient
df$pt_cluster = paste(df$PATIENT_ID, "_", df$cluster_id, sep = "")
df = unique(select(df, pt_cluster, symbol, cellular_prevalence))
df$cellular_prevalence = as.numeric(df$cellular_prevalence)

# since there are some cases where there are multiple mutations in the same gene in the same cluster, we need to annotate the symbol based on the order to be able to use dcast

unique_clusters = as.data.frame(unique(df$pt_cluster))

for(i in 1:nrow(unique_clusters)){
  df_sub = subset(df, df$pt_cluster == unique_clusters$`unique(df$pt_cluster)`[i])
  
  df_dup = df_sub[duplicated(df_sub$symbol),]
  
  if(nrow(df_dup) > 0){
    gene_dup = df_dup$symbol[1]
    
    n = 1
    
    for(j in 1:nrow(df)){
      
      if(df$pt_cluster[j] == unique_clusters$`unique(df$pt_cluster)`[i] & df$symbol[j] == gene_dup){
        
        df$symbol[j] = paste(gene_dup, "_", n, sep = "")
        
        n = n + 1
        
      }
      if(n == 2){
        if(df$pt_cluster[j] == unique_clusters$`unique(df$pt_cluster)`[i] & df$symbol[j] == gene_dup){
          
          k = 2
          
          df$symbol[j] = paste(gene_dup, "_", k, sep = "")
          
        } 
      }
    }
  }
}

# refomate the data for use in clustering heatmap

# plot raw PyClone results (gene prevelance by clone)
df_2 = reshape2::dcast(df, pt_cluster ~ symbol, value.var="cellular_prevalence")

rownames(df_2) <- df_2$pt_cluster
df_2$pt_cluster <- NULL
df_2[is.na(df_2)] <- 0

df_2 = as.data.frame(t(data.matrix(df_2)))

# cluster and viauslize the data
p = pheatmap(df_2, 
             color = colorRampPalette(c( "#2d004b", "#FDE257FF"))(25),
             legend_labels = c("0 %", "25 %", "50 %", "75 %", "100 %", "Cellular\nPrevalence"),
             border_color = NA,
             show_colnames = F,
             cluster_rows = T,
             cluster_cols = T,
             fontsize = 5,
             treeheight_row = 10,
             treeheight_col = 25
)
ggsave(p, filename = "~/Desktop/MetaAML_results/Figure_5/cluster_by_clones.png", dpi = 300, width = 4, height = 3.5, units = "in")



# plot per patient PyClone results (gene prevelance by patient)
df2 = df

df2$pt_cluster = gsub("_.*","",df2$pt_cluster)

pts = as.data.frame(unique(df2$pt_cluster))

pt_gene = list()
a = 1

for(i in 1:nrow(pts)){
  print(i)
  pt_sub = subset(df2, pt_cluster == pts$`unique(df2$pt_cluster)`[i])
  genes = as.data.frame(unique(pt_sub$symbol))
  for(j in 1:nrow(genes)){
    print(j)
    pt_gene_sub = subset(pt_sub, symbol == genes$`unique(pt_sub$symbol)`[j])
    if(nrow(pt_gene_sub) == 1){
      pt_gene[[a]] = pt_gene_sub 
      a = a + 1
    }
    if(nrow(pt_gene_sub) == 2){
      pt_gene_sub = pt_gene_sub %>% 
        arrange(desc(cellular_prevalence))
      
      pt_gene_sub$symbol[2] = paste(genes$`unique(pt_sub$symbol)`[j], "_2", sep = "")

      pt_gene[[a]] = pt_gene_sub 
      a = a + 1
    }
    if(nrow(pt_gene_sub) == 3){
      pt_gene_sub = pt_gene_sub %>% 
        arrange(desc(cellular_prevalence))
      
      pt_gene_sub$symbol[2] = paste(genes$`unique(pt_sub$symbol)`[j], "_2", sep = "")
      pt_gene_sub$symbol[3] = paste(genes$`unique(pt_sub$symbol)`[j], "_3", sep = "")
      
      pt_gene[[a]] = pt_gene_sub 
      a = a + 1
    }
  }
}
pt_gene_list <- do.call(rbind, pt_gene)

# test to make sure there are no duplicate patient-gene pairs in order for dcast to work
for (i in 1:nrow(pt_gene_list)) {
  pt_gene_list$not_unique[i] <- length(which(pt_gene_list$pt_cluster == pt_gene_list$pt_cluster[i] & pt_gene_list$symbol == pt_gene_list$symbol[i]))
}

# shange to wide format
df2 = reshape2::dcast(pt_gene_list, pt_cluster ~ symbol, value.var="cellular_prevalence")

rownames(df2) <- df2$pt_cluster
df2$pt_cluster <- NULL
df2[is.na(df2)] <- 0

df2 = t(data.matrix(df2))

# plot
# add annotations
# Data frame with column annotations
df = read_tsv("~/Desktop/MetaAML_results/Figure_5/PyClone_final.tsv")

df_annotations = select(df, PATIENT_ID, Risk, isDenovo, isRelapse, isTransformed, Cohort)

df_annotations$Subset = NA

for(i in 1:nrow(df_annotations)){
  if (df_annotations$isDenovo[i] == "TRUE" & df_annotations$isRelapse[i] == "FALSE") {
    df_annotations$Subset[i] <- "De novo"
  }
  if (df_annotations$isDenovo[i] == "TRUE" & df_annotations$isRelapse[i] == "TRUE") {
    df_annotations$Subset[i] <- "Relapse"
  }
  if (df_annotations$isRelapse[i] == "TRUE" & df_annotations$isTransformed[i] == "FALSE") {
    df_annotations$Subset[i] <- "Relapse"
  }
  if (df_annotations$isRelapse[i] == "TRUE" & df_annotations$isTransformed[i] == "TRUE") {
    df_annotations$Subset[i] <- "Relapse"
  }
  if (df_annotations$isTransformed[i] == "TRUE" & df_annotations$isRelapse[i] == "FALSE") {
    df_annotations$Subset[i] <- "Transformed"
  }
  if (df_annotations$isDenovo[i] == "FALSE" & df_annotations$isRelapse[i] == "FALSE" & df_annotations$isTransformed[i] == "FALSE") {
    df_annotations$Subset[i] <- "Other"
  }
}
for(i in 1:nrow(df_annotations)){
  if (df_annotations$Risk[i] == "FavorableOrIntermediate" | df_annotations$Risk[i] == "IntermediateOrAdverse"){
    df_annotations$Risk[i] = "Intermediate"
  }
  if (df_annotations$Risk[i] == "Good"){
    df_annotations$Risk[i] = "Favorable"
  }
  if (df_annotations$Risk[i] == "Poor"){
    df_annotations$Risk[i] = "Adverse"
  }
  if (df_annotations$Risk[i] == "N.D."){
    df_annotations$Risk[i] = "Unknown"
  }
}

for(i in 1:nrow(df_annotations)){
  if(df_annotations$Cohort[i] == "TCGA"){
    df_annotations$Subset[i] = "De novo"
  }
}

df_annotations = as.data.frame(unique(na.omit(select(df_annotations, PATIENT_ID, Risk, Subset, Cohort))))

# mat_col <- data.frame(group = col_groups)
rownames(df_annotations) <- df_annotations$PATIENT_ID
df_annotations$PATIENT_ID = NULL

# List with colors for each annotation.
col = list(Risk = c("Adverse" = "#E64B35FF",
                    "Intermediate" = "#8491B4FF",
                    "Favorable" = "#00A087FF", 
                    "Unknown" = "#767676FF"),
           Subset = c("De novo" = "#C16622FF", 
                      "Transformed" = "#767676FF", 
                      "Relapse" = "#800000FF", 
                      "Other" = "#FFA319FF"),
           Cohort = c("Tyner" = "#0073C2FF", 
                      "TCGA" = '#EFC000FF'))

p = pheatmap(df2, 
             color = colorRampPalette(c("white", "darkred"))(50),
             annotation_col = df_annotations,
             annotation_colors = col,
             legend_breaks = c(0, 0.25, 0.5, 0.75, 1, max(df2)), 
             main = "", legend_labels = c("0 %", "25 %", "50 %", "75 %", "100 %", "Cellular\nPrevalence"),
             legend = TRUE,
             border_color = NA,
             show_colnames = F,
             cluster_rows = T,
             cluster_cols = T,
             fontsize = 5,
             treeheight_row = 5,
             treeheight_col = 10
)
ggsave(p, filename = "~/Desktop/MetaAML_results/Figure_5/cluster_pts_by_gene_celular_prevalence.png", dpi = 300, width = 6, height = 5, units = "in")

# for visualization purposes, make heatmaps with and without the second mutations in the same genes. For simple heatmap, only include the mutation with the highest cellular prevalence
df2_simple = subset(pt_gene_list, !grepl("_", pt_gene_list$symbol))
# shange to wide format
df2 = reshape2::dcast(df2_simple, pt_cluster ~ symbol, value.var="cellular_prevalence")

rownames(df2) <- df2$pt_cluster
df2$pt_cluster <- NULL
df2[is.na(df2)] <- 0

df2 = t(data.matrix(df2))

# plot
# add annotations
# Data frame with column annotations
df = read_tsv("~/Desktop/MetaAML_results/Figure_5/PyClone_final.tsv")

df_annotations = select(df, PATIENT_ID, Risk, isDenovo, isRelapse, isTransformed, Cohort)

df_annotations$Subset = NA

for(i in 1:nrow(df_annotations)){
  if (df_annotations$isDenovo[i] == "TRUE" & df_annotations$isRelapse[i] == "FALSE") {
    df_annotations$Subset[i] <- "De novo"
  }
  if (df_annotations$isDenovo[i] == "TRUE" & df_annotations$isRelapse[i] == "TRUE") {
    df_annotations$Subset[i] <- "Relapse"
  }
  if (df_annotations$isRelapse[i] == "TRUE" & df_annotations$isTransformed[i] == "FALSE") {
    df_annotations$Subset[i] <- "Relapse"
  }
  if (df_annotations$isRelapse[i] == "TRUE" & df_annotations$isTransformed[i] == "TRUE") {
    df_annotations$Subset[i] <- "Relapse"
  }
  if (df_annotations$isTransformed[i] == "TRUE" & df_annotations$isRelapse[i] == "FALSE") {
    df_annotations$Subset[i] <- "Transformed"
  }
  if (df_annotations$isDenovo[i] == "FALSE" & df_annotations$isRelapse[i] == "FALSE" & df_annotations$isTransformed[i] == "FALSE") {
    df_annotations$Subset[i] <- "Other"
  }
}
for(i in 1:nrow(df_annotations)){
  if (df_annotations$Risk[i] == "FavorableOrIntermediate" | df_annotations$Risk[i] == "IntermediateOrAdverse"){
    df_annotations$Risk[i] = "Intermediate"
  }
  if (df_annotations$Risk[i] == "Good"){
    df_annotations$Risk[i] = "Favorable"
  }
  if (df_annotations$Risk[i] == "Poor"){
    df_annotations$Risk[i] = "Adverse"
  }
  if (df_annotations$Risk[i] == "N.D."){
    df_annotations$Risk[i] = "Unknown"
  }
}

for(i in 1:nrow(df_annotations)){
  if(df_annotations$Cohort[i] == "TCGA"){
    df_annotations$Subset[i] = "De novo"
  }
}

df_annotations = as.data.frame(unique(na.omit(select(df_annotations, PATIENT_ID, Risk, Subset, Cohort))))

# mat_col <- data.frame(group = col_groups)
rownames(df_annotations) <- df_annotations$PATIENT_ID
df_annotations$PATIENT_ID = NULL

# List with colors for each annotation.
col = list(Risk = c("Adverse" = "#E64B35FF",
                    "Intermediate" = "#8491B4FF",
                    "Favorable" = "#00A087FF", 
                    "Unknown" = "#767676FF"),
           Subset = c("De novo" = "#C16622FF", 
                      "Transformed" = "#767676FF", 
                      "Relapse" = "#800000FF", 
                      "Other" = "#FFA319FF"),
           Cohort = c("Tyner" = "#0073C2FF", 
                      "TCGA" = '#EFC000FF'))

p = pheatmap(df2, 
             color = colorRampPalette(c("white", "darkred"))(50),
             annotation_col = df_annotations,
             annotation_colors = col,
             legend_breaks = c(0, 0.25, 0.5, 0.75, 1, max(df2)), 
             main = "", legend_labels = c("0 %", "25 %", "50 %", "75 %", "100 %", "Cellular\nPrevalence"),
             legend = TRUE,
             border_color = NA,
             show_colnames = F,
             cluster_rows = T,
             cluster_cols = T,
             fontsize = 5,
             treeheight_row = 5,
             treeheight_col = 10 
)
ggsave(p, filename = "~/Desktop/MetaAML_results/Figure_5/cluster_pts_by_gene_celular_prevalence_no_bi_allelic.png", dpi = 300, width = 4
       , height = 3, units = "in")

 # A. Histogram of number of clones
# B. Histogram of mutation burden
# C. Histogram of ages
# D. Histogram of SDI
# E. Mutation burden vs. clonality
# F. Clonality vs. age
# G. Mutation burden vs. age
# H. Overall Survival
# I. Survival, two clone bins
# J. Survival, four clone bins
# K. Survival, Shannon and mutation for 15 bins
# L. Survival, Shannon and mutation analysis
# M. Shannon scatterplots
# N. Range of VAFs within each clone
# O. Math Analysis
# P. Clonal Evolution Analysis

dir.create("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs")
setwd("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs")
####################################
# A. Histogram of number of clones #
####################################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones, total_mutations) %>%
  distinct() %>%
  filter(!is.na(number_of_clones)) 
# 358 samples

data_for_plot %>%
  ggplot(aes(number_of_clones)) + geom_histogram(binwidth = 1, fill = "lightblue", color = "black") + ggtitle("Frequency of Clonal Populations") + xlab("Number of Clones") + ylab("Frequency") + theme_cowplot() + ggsave("Clonality_histogram.pdf", units="in", width=5, height=4, dpi=300)

###################################
# B. Histogram of mutation burden #
###################################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations) %>%
  distinct() %>%
  filter(!is.na(total_mutations))
# 358 samples

data_for_plot %>%
  ggplot(aes(total_mutations)) + geom_histogram(binwidth = 1, fill = "lightblue", color = "black") + ggtitle("Frequency of Mutation Burdens") + xlab("Number of Mutations") + ylab("Frequency") + theme_cowplot() + ggsave("Mutation_burden_histogram.pdf", units="in", width=5, height=4, dpi=300)

########################
# C. Histogram of ages #
########################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE) %>%
  distinct() %>%
  filter(!is.na(AGE))
# 358 samples

data_for_plot %>%
  ggplot(aes(AGE)) + geom_histogram(binwidth = 1, fill = "lightblue", color = "black") + ggtitle("Frequency of Ages") + xlab("Age") + ylab("Frequency") + theme_cowplot() + ggsave("Age_histogram.pdf", units="in", width=5, height=4, dpi=300)

####################################
# D. Histogram of SDI              #
####################################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, SDI) %>%
  distinct() %>%
  filter(!is.na(SDI))
# 313 samples

data_for_plot %>%
  ggplot(aes(SDI)) + geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") + ggtitle("Frequency of Shannon Diversity Indeces") + xlab("SDI") + ylab("Frequency") + theme_cowplot() + ggsave("SDI_histogram.pdf", units="in", width=5, height=4, dpi=300)

####################################
# E. Mutation burden vs. clonality #
####################################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones, total_mutations) %>%
  distinct() %>%
  filter(!is.na(number_of_clones)) %>%
  filter(!is.na(total_mutations))
# 358 samples

data_for_plot %>%
  ggplot(aes(number_of_clones, total_mutations)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mutation Burden vs. Clonality") + xlab("Number of Clones") + ylab("Number of Mutations") + stat_regline_equation(label.x = 10) + stat_cor(label.x = 10, label.y = 50) + theme_cowplot() + ggsave("Mutation_burden_vs_clonality.pdf", units="in", width=5, height=4, dpi=300)

########################
# F. Clonality vs. Age #
########################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones, AGE) %>%
  distinct() %>%
  filter(!is.na(number_of_clones)) %>%
  filter(!is.na(AGE))
# 358 samples

data_for_plot %>%
  ggplot(aes(AGE, number_of_clones)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Clonality vs. Age") + xlab("Age") + ylab("Number of Clones") + stat_regline_equation() + stat_cor(label.y = 17) + theme_cowplot() + ggsave("Clonality_vs_age.pdf", units="in", width=5, height=4, dpi=300)

##############################
# G. Mutation Burden vs. Age #
##############################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations, AGE) %>%
  distinct() %>%
  filter(!is.na(total_mutations)) %>%
  filter(!is.na(AGE))
# 358 samples

data_for_plot %>%
  ggplot(aes(AGE, total_mutations)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mutation Burden vs. Age") + xlab("Age") + ylab("Number of Mutations") + stat_regline_equation(label.x = 65) + stat_cor(label.x = 65, label.y = 55) + theme_cowplot() + ggsave("Mutation_burden_vs_age.pdf", units="in", width=5, height=4, dpi=300)

#######################
# H. Overall Survival #
#######################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones, SDI, total_mutations) %>%
  distinct() %>%
  filter(!is.na(Censor)) %>%
  filter(!is.na(overallSurvival)) %>%
  filter(!is.na(number_of_clones)) %>%
  filter(!is.na(total_mutations))
# 295 samples

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$one_bin <- 1

fit <- survfit(Surv(overallSurvival, Censor) ~ one_bin, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = "black") + ggsave("Overall_survival.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Overall Survival")


## add in mutation burden to analysis
mean_mutation <- mean(data_for_plot$total_mutations)
data_for_plot$mutation_bin <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] < mean_mutation) {
    data_for_plot$mutation_bin[i] <- "Low"
  } else {
    data_for_plot$mutation_bin[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ mutation_bin, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Overall_survival_with_mutation.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Overall Survival with Mutation High/low Bins")

###############################
# I. Survival, two clone bins #
###############################
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= 6) {
    data_for_plot$two_bins[i] <- "1-6"
  } else {
    data_for_plot$two_bins[i] <- "7+"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_two_bins.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Survival with 2 Clone Bins")

## add in mutation burden to analysis
data_for_plot <- data_for_plot %>%
  group_by(two_bins) %>%
  mutate(mean_mutation_two_bins = mean(total_mutations)) %>%
  ungroup() %>%
  mutate(two_bins = fct_relevel(two_bins, "1-10", "11+"))

# assign high or low mutation by bin
data_for_plot$mutation_bin_two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] < data_for_plot$mean_mutation_two_bins[i]) {
    data_for_plot$mutation_bin_two_bins[i] <- "Low"
  } else {
    data_for_plot$mutation_bin_two_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ mutation_bin_two_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_two_bins_with_mutation.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Mutation High/low within 2 Clone Bins")

## boxplot of mutation burden
data_for_plot %>%
  ggplot(aes(two_bins, total_mutations)) + geom_boxplot() + xlab("Number of Clones") + ggtitle("Mutation Burden by Clone Bin") + ylab("Total Mutation Burden") + theme_cowplot() + ggsave("Mutation_boxplot_two_clone_bins.pdf", units="in", width=5, height=4, dpi=300)

################################
# J. Survival, four clone bins #
################################
data_for_plot$four_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= 3) {
    data_for_plot$four_bins[i] <- "1-3"
  } else if (data_for_plot$number_of_clones[i] <= 6) {
    data_for_plot$four_bins[i] <- "4-6"
  } else if (data_for_plot$number_of_clones[i] <= 9) {
    data_for_plot$four_bins[i] <- "7-9"
  } else {
    data_for_plot$four_bins[i] <- "10+"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ four_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue", "red", "gold")) + ggsave("Survival_four_bins.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Survival with 4 Clone Bins")

## add in mutation burden to analysis
data_for_plot <- data_for_plot %>%
  group_by(four_bins) %>%
  mutate(mean_mutation_four_bins = mean(total_mutations)) %>%
  ungroup() %>%
  mutate(four_bins = fct_relevel(four_bins, "1-5", "6-10", "11-15", "16+"))

# assign high or low mutation burden by bin
data_for_plot$mutation_bin_four_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] < data_for_plot$mean_mutation_four_bins[i]) {
    data_for_plot$mutation_bin_four_bins[i] <- "Low"
  } else {
    data_for_plot$mutation_bin_four_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ mutation_bin_four_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_four_bins_with_mutation.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Mutation High/low within 4 Clone Bins")

## boxplot of mutation burden
data_for_plot %>%
  ggplot(aes(four_bins, total_mutations)) + geom_boxplot() + xlab("Number of Clones") + ggtitle("Mutation Burden by Clone Bin") + ylab("Total Mutation Burden") + theme_cowplot() + ggsave("Mutation_boxplot_four_clone_bins.pdf", units="in", width=5, height=4, dpi=300)

## Mutation high and low for 10 clone bins
data_for_plot$ten_bins <- NA

for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= 9) {
    data_for_plot$ten_bins[i] <- as.character(data_for_plot$number_of_clones[i])
  } else {
    data_for_plot$ten_bins[i] <- "10-13"
  }
}

## Add mutation burden to analysis
data_for_plot <- data_for_plot %>%
  group_by(ten_bins) %>%
  mutate(mean_mutation_ten_bins = mean(total_mutations)) %>%
  ungroup() %>%
  mutate(ten_bins = fct_relevel(ten_bins, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10-13"))

# assign high or low mutation burden by bin
data_for_plot$mutation_bin_ten_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] < data_for_plot$mean_mutation_ten_bins[i]) {
    data_for_plot$mutation_bin_ten_bins[i] <- "Low"
  } else {
    data_for_plot$mutation_bin_ten_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ mutation_bin_ten_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_ten_bins_with_mutation.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Mutation High/low within 10 Clone Bins")

## boxplot of mutation burden for ten clone bins
data_for_plot %>%
  ggplot(aes(ten_bins, total_mutations)) + geom_boxplot() + xlab("Number of Clones") + ggtitle("Mutation Burden by Clone Bin") + ylab("Mutation Burden") + theme_cowplot() + ggsave("Mutation_boxplot_ten_clone_bins.pdf", units="in", width=5, height=4, dpi=300)

######################################
# K. Shannon high/low for clone bins #
######################################
data_for_plot <- data_for_plot %>%
  filter(!is.na(SDI))

## Shannon high/low for overall survival
mean_shannon <- mean(data_for_plot$SDI)
data_for_plot$shannon_bin <- NA
for (i in 1:nrow(data_for_plot)) {
  if (is.na(data_for_plot$SDI[i])) {
    data_for_plot$shannon_bin[i] <- NA
  } else if (data_for_plot$SDI[i] < mean_shannon) {
    data_for_plot$shannon_bin[i] <- "Low"
  } else {
    data_for_plot$shannon_bin[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ shannon_bin, data = data_for_plot)
print(fit)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata",  palette = c(cbPalette[5], cbPalette[6]), legend.labs = c("Shannon High", "Shannon Low"), legend.title = "", legend = c(0.4, 0.9)) + ggsave("Overall_survival_with_shannon.pdf", units="in", width=4, height=3, dpi=300)

## Shannon high/low for two clone bins
data_for_plot <- data_for_plot %>%
  group_by(two_bins) %>%
  mutate(mean_shannon_two_bins = mean(SDI))

# assign high or low SDI by bin
data_for_plot$shannon_bin_two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$SDI[i] < data_for_plot$mean_shannon_two_bins[i]) {
    data_for_plot$shannon_bin_two_bins[i] <- "Low"
  } else {
    data_for_plot$shannon_bin_two_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ shannon_bin_two_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years",  risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_two_bins_with_shannon.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Shannon High/low within 2 Clone Bins")

## boxplot of SDI for two clone bins
data_for_plot %>%
  ggplot(aes(two_bins, SDI)) + geom_boxplot() + xlab("Number of Clones") + ggtitle("Shannon Diversity Index by Clone Bin") + ylab("Shannon Diversity Index") + theme_cowplot() + ggsave("SDI_boxplot_four_clone_bins.pdf", units="in", width=5, height=4, dpi=300)

## Shannon high/low for four bins
data_for_plot <- data_for_plot %>%
  group_by(four_bins) %>%
  mutate(mean_shannon_four_bins = mean(SDI))

# assign high or low SDI by bin
data_for_plot$shannon_bin_four_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$SDI[i] < data_for_plot$mean_shannon_four_bins[i]) {
    data_for_plot$shannon_bin_four_bins[i] <- "Low"
  } else {
    data_for_plot$shannon_bin_four_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ shannon_bin_four_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years",  risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_four_bins_with_shannon.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Shannon High/low within 4 Clone Bins")

## boxplot of SDI for four clone bins
data_for_plot %>%
  ggplot(aes(four_bins, SDI)) + geom_boxplot() + xlab("Number of Clones") + ggtitle("Shannon Diversity Index by Clone Bin") + ylab("Shannon Diversity Index") + theme_cowplot() + ggsave("SDI_boxplot_four_clone_bins.pdf", units="in", width=5, height=4, dpi=300)


## Shannon high/low for 10 bins
data_for_plot <- data_for_plot %>%
  group_by(ten_bins) %>%
  mutate(mean_shannon_ten_bins = mean(SDI)) %>%
  ungroup()

# assign high or low SDI by bin
data_for_plot$shannon_bin_ten_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$SDI[i] < data_for_plot$mean_shannon_ten_bins[i]) {
    data_for_plot$shannon_bin_ten_bins[i] <- "Low"
  } else {
    data_for_plot$shannon_bin_ten_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ shannon_bin_ten_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_ten_bins_with_shannon.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Shannon High/low within 10 Clone Bins")

## boxplot of SDI for ten clone bins
data_for_plot %>%
  ggplot(aes(ten_bins, SDI)) + geom_boxplot() + xlab("Number of Clones") + ggtitle("Shannon Diversity Index by Clone Bin") + ylab("Shannon Diversity Index") + theme_cowplot() + ggsave("SDI_boxplot_ten_clone_bins.pdf", units="in", width=5, height=4, dpi=300)




####################################
# L. Shannon and mutation analysis #
####################################
## shannon high and low by mutation bin
data_for_plot$four_mutation_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] <= 10) {
    data_for_plot$four_mutation_bins[i] <- "1-10"
  } else if (data_for_plot$total_mutations[i] <= 20) {
    data_for_plot$four_mutation_bins[i] <- "11-20"
  } else if (data_for_plot$total_mutations[i] <= 30) {
    data_for_plot$four_mutation_bins[i] <- "21-30"
  } else {
    data_for_plot$four_mutation_bins[i] <- "31+"
  }
}

# Add Shannon to analysis
data_for_plot <- data_for_plot %>%
  group_by(four_mutation_bins) %>%
  mutate(mean_shannon_four_mutation_bins = mean(SDI))

# assign high or low Shannon by bin
data_for_plot$shannon_bin_four_mutation_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$SDI[i] < data_for_plot$mean_shannon_four_mutation_bins[i]) {
    data_for_plot$shannon_bin_four_mutation_bins[i] <- "Low"
  } else {
    data_for_plot$shannon_bin_four_mutation_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ shannon_bin_four_mutation_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_mutation_bins_with_shannon.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Shannon High/low within 4 Mutation Bins")

## boxplot of SDI for four mutation bins
data_for_plot %>%
  ggplot(aes(four_mutation_bins, SDI)) + geom_boxplot() + xlab("Mutation Bin") + ggtitle("Shannon Diversity Index by Mutation Bin") + ylab("Shannon Diversity Index") + theme_cowplot() + ggsave("SDI_boxplot_four_mutation_bins.pdf", units="in", width=5, height=4, dpi=300)


## mutation high and low by shannon bin (five bins)
data_for_plot$five_shannon_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$SDI[i] <= 1) {
    data_for_plot$five_shannon_bins[i] <- "0.0-1.0"
  } else if (data_for_plot$SDI[i] <= 1.5) {
    data_for_plot$five_shannon_bins[i] <- "1.0-1.5"
  } else if (data_for_plot$SDI[i] <= 2) {
    data_for_plot$five_shannon_bins[i] <- "1.5-2.0"
  } else if (data_for_plot$SDI[i] <= 2.5) {
    data_for_plot$five_shannon_bins[i] <- "2.0-2.5"
  } else {
    data_for_plot$five_shannon_bins[i] <- "2.5-3.0"
  }
}

# Add mutation burden to analysis
data_for_plot <- data_for_plot %>%
  group_by(five_shannon_bins) %>%
  mutate(mean_mutation_five_shannon_bins = mean(total_mutations))

# assign high or low Shannon by bin
data_for_plot$mutation_bin_five_shannon_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] < data_for_plot$mean_mutation_five_shannon_bins[i]) {
    data_for_plot$mutation_bin_five_shannon_bins[i] <- "Low"
  } else {
    data_for_plot$mutation_bin_five_shannon_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ mutation_bin_five_shannon_bins, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Survival_shannon_bins_with_mutation.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Mutation High/low within 5 Shannon Bins")

## boxplot of mutation burden for five SDI bins
data_for_plot %>%
  ggplot(aes(five_shannon_bins, SDI)) + geom_boxplot() + xlab("SDI Bin") + ggtitle("Mutation Burden by SDI bin") + ylab("Mutation Burden") + theme_cowplot() + ggsave("Mutation_boxplot_five_SDI_bins.pdf", units="in", width=5, height=4, dpi=300)

###########################
# M. Shannon scatterplots #
###########################

## Shannon vs. age
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, SDI) %>%
  distinct() %>%
  filter(!is.na(AGE)) %>%
  filter(!is.na(SDI))

data_for_plot %>%
  ggplot(aes(AGE, SDI)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Shannon Diversity Index vs. Age") + xlab("Age") + ylab("Shannon Diversity Index") + stat_regline_equation() + stat_cor(label.y = 2.55) + theme_cowplot() + ggsave("SDI_vs_Age.pdf", units = "in", width=5, height=4, dpi=300)

## Shannon vs. mutation burden
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations, SDI) %>%
  distinct() %>%
  filter(!is.na(total_mutations)) %>%
  filter(!is.na(SDI))

data_for_plot %>%
  ggplot(aes(total_mutations, SDI)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Shannon Diversity Index vs. Mutation Burden") + xlab("Mutation Burden") + ylab("Shannon Diversity Index") + stat_regline_equation(label.y = 3.75) + stat_cor(label.y = 3.25)+ theme_cowplot() + ggsave("SDI_vs_mutation.pdf", units = "in", width=5, height=4, dpi=300)


## Shannon vs. clonality
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones, SDI) %>%
  distinct() %>%
  filter(!is.na(number_of_clones)) %>%
  filter(!is.na(SDI))

data_for_plot %>%
  ggplot(aes(number_of_clones, SDI)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Shannon Diversity Index vs. Clonality") + xlab("Number of Clones") + ylab("Shannon Diversity Index") + stat_regline_equation(label.y = 3.75) + stat_cor(label.y = 3.25) + theme_cowplot() + ggsave("SDI_vs_clonality.pdf", units = "in", width=5, height=4, dpi=300)

######################################
# N. Range of VAFs within each clone #
######################################
## For clones that have more than one mutation, create a histogram of the max-min VAF within each clone to get an idea of the spread of VAF sizes that can still be clustered into the same clone by PyClone
data_for_plot <- PyClone %>%
  select(PATIENT_ID, cluster_id, variant_allele_frequency) %>%
  mutate(id = str_c(PATIENT_ID, as.character(cluster_id), sep = "_")) %>%
  group_by(id) %>%
  summarise(range = max(variant_allele_frequency) - min(variant_allele_frequency)) %>%
  filter(range > 0) %>%
  distinct()

data_for_plot %>%
  ggplot(aes(range)) + geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") + ggtitle("Ranges in VAFs within each clone of size >= 2") + xlab("VAF Range") + ylab("Frequency") + theme_cowplot() + ggsave("VAF_range_histogram.pdf", units = "in", width=5, height=4, dpi=300)

## Repeat this analysis but with the standard deviation of the VAFs within each clone of >=2
data_for_plot <- PyClone %>%
  select(PATIENT_ID, cluster_id, variant_allele_frequency) %>%
  mutate(id = str_c(PATIENT_ID, as.character(cluster_id), sep = "_")) %>%
  group_by(id) %>%
  summarise(std = sd(variant_allele_frequency)) %>%
  filter(std > 0) %>%
  distinct()

data_for_plot %>%
  ggplot(aes(std)) + geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") + ggtitle("VAF standard deviations within each clone of size >= 2") + xlab("VAF Standard Deviation") + ylab("Frequency") + theme_cowplot() + ggsave("VAF_std_histogram.pdf", units = "in", width=5, height=4, dpi=300)

## If the number of clones is greater than one, for each entry, determine the minimum difference between the VAF of that entry and any VAF in a different clone to determine how granular PyClone is
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones, variant_allele_frequency, cluster_id) %>%
  filter(number_of_clones > 1) %>%
  filter(!is.na(variant_allele_frequency)) %>%
  filter(!is.na(cluster_id))

data_for_plot$min_diff <- NA

for (i in 1:nrow(data_for_plot)) {
  data_for_plot$min_diff[i] <- abs(min(data_for_plot$variant_allele_frequency[which(data_for_plot$PATIENT_ID == data_for_plot$PATIENT_ID[i] & data_for_plot$cluster_id != data_for_plot$cluster_id[i])] - data_for_plot$variant_allele_frequency[i]))
}

data_for_plot %>%
  filter(min_diff > 0 & min_diff != Inf) %>%
  ggplot(aes(min_diff)) + geom_histogram(color = "black", fill = "lightblue", binwidth = 0.01) + ggtitle("Minimum difference between VAFs of different clones") + xlab("Minimum VAF Difference") + ylab("Frequency") + theme_cowplot() + ggsave("Min_VAF_Difference_Histogram.pdf", units = "in", width = 5, height = 4, dpi = 300)

####################
# O. Math Analysis #
####################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations, number_of_clones, variant_allele_frequency, overallSurvival, Censor, SDI, AGE) %>%
  filter(total_mutations > 1) %>%
  filter(!is.na(Censor)) %>%
  filter(!is.na(overallSurvival)) %>%
  group_by(PATIENT_ID) %>%
  mutate(MATH_Score = 100*mad(variant_allele_frequency)/median(variant_allele_frequency)) %>%
  select(PATIENT_ID, total_mutations, number_of_clones, overallSurvival, Censor, MATH_Score, SDI, AGE) %>%
  filter(complete.cases(MATH_Score)) %>%
  distinct()
# 269 observations

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)

## histogram of MATH scores
data_for_plot %>%
  ggplot(aes(MATH_Score)) + geom_histogram(binwidth = 1, fill = "lightblue", color = "black") + ggtitle("Frequency of MATH Scores") + xlab("MATH Score") + ylab("Frequency") + theme_cowplot() + ggsave("MATH_histogram.pdf", units="in", width=5, height=4, dpi=300)

## Scatterplot of MATH vs. Number of mutations
data_for_plot %>%
  ggplot(aes(total_mutations, MATH_Score)) + geom_point() + geom_smooth(method = "lm") + ggtitle("MATH Score vs. Mutation Burden") + xlab("Number of Mutations") + ylab("MATH Score") + stat_regline_equation(label.x = 15) + stat_cor(label.x = 15, label.y = 80) + theme_cowplot() + ggsave("MATH_vs_mutation_burden.pdf", units="in", width=5, height=4, dpi=300)

## Scatterplot of MATH vs. number of clones
data_for_plot %>%
  ggplot(aes(number_of_clones, MATH_Score)) + geom_point() + geom_smooth(method = "lm") + ggtitle("MATH Score vs. Clonality") + xlab("Number of Clones") + ylab("MATH Score") + stat_regline_equation(label.x = 5) + stat_cor(label.x = 5, label.y = 80) + theme_cowplot() + ggsave("MATH_vs_clonality.pdf", units="in", width=5, height=4, dpi=300)

## Scatterplot of MATH vs. Age
data_for_plot %>%
  ggplot(aes(AGE, MATH_Score)) + geom_point() + geom_smooth(method = "lm") + ggtitle("MATH Score vs. AGE") + xlab("Age") + ylab("MATH Score") + stat_regline_equation(label.x = 65) + stat_cor(label.x = 65, label.y = 80) + theme_cowplot() + ggsave("MATH_vs_age.pdf", units="in", width=5, height=4, dpi=300)

## Scatterplot of MATH vs. SDI
data_for_plot %>%
  ggplot(aes(SDI, MATH_Score)) + geom_point() + geom_smooth(method = "lm") + ggtitle("MATH Score vs. SDI") + xlab("SDI") + ylab("MATH Score") + stat_regline_equation(label.x = 1) + stat_cor(label.x = 1, label.y = 80) + theme_cowplot() + ggsave("MATH_vs_SDI.pdf", units="in", width=5, height=4, dpi=300)

## Survival with MATH
mean_math <- mean(data_for_plot$MATH_Score, na.rm = TRUE)
data_for_plot <- data_for_plot %>%
  mutate(MATH_bin = ifelse(MATH_Score < mean_math, "low", "high"))

fit <- survfit(Surv(overallSurvival, Censor) ~ MATH_bin, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", palette = "jco", legend.labs = c("MATH High", "MATH Low"), legend.title = "", legend = c(0.4, 0.9)) + ggsave("Overall_Survival_with_MATH_bins.pdf", units="in", width=4, height=3, dpi=300)


# ClonEvol ####
################################
# P. Clonal Evolution Analysis #
################################
data_for_plot <- PyClone_with_all_mutations %>%
  select(PATIENT_ID, cluster_id, cellular_prevalence) %>%
  na.omit() %>%
  mutate(cellular_prevalence = cellular_prevalence*100)

data_for_plot$cluster_id = data_for_plot$cluster_id+1

dir.create("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Evolution/")

for (i in 1:length(unique(data_for_plot$PATIENT_ID))){
  # for (i in 262:262) {
  patient <- unique(data_for_plot$PATIENT_ID)[i]
  dir.create(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Evolution/", patient))
  setwd(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Evolution/", patient))
  height <- length(which(data_for_plot$PATIENT_ID == patient))
  x <- data.frame(cluster = rep(NA, height), bogus.vaf = 0, second.vaf = 0)
  input_data <- filter(data_for_plot, PATIENT_ID == patient)
  
  for (i in 1:height) {
    x$cluster[i] <- input_data$cluster_id[i]
    x[i, 2] <- input_data$cellular_prevalence[i]
    x[i, 3] <- input_data$cellular_prevalence[i]
  }
  vaf.col.names <- grep('.vaf', colnames(x), value=T)
  sample.names <- gsub('.vaf', '', vaf.col.names)
  x[, sample.names] <- x[, vaf.col.names]
  vaf.col.names <- sample.names
  sample.groups <- c("bogus", "second")
  names(sample.groups) <- vaf.col.names
  x <- x[order(x$cluster),]
  clone.colors <- NULL
  pdf('box.pdf', width = 3, height = 3, useDingbats = FALSE, title = '')
  plot.variant.clusters(x, cluster.col.name = 'cluster', show.cluster.size = TRUE, cluster.size.text.color = 'blue', vaf.col.names = vaf.col.names, vaf.limits = 100, sample.title.size = 20, box = FALSE, jitter = TRUE, jitter.shape = 1, jitter.color = "black", jitter.center.method = 'median', jitter.alpha = 1)
  dev.off()
  
  x<- x %>%
    group_by(cluster) %>%
    mutate(average_vaf = mean(bogus))
  founding.cluster <- x$cluster[first(which(x$average_vaf == max(x$average_vaf)))]
  
  y <- infer.clonal.models(variants = x, cluster.col.name = 'cluster', vaf.col.names = vaf.col.names, sample.groups = sample.groups, cancer.initiation.model = 'monoclonal', subclonal.test = 'bootstrap',  subclonal.test.model = 'non-parametric', num.boots = 1000, founding.cluster = founding.cluster, cluster.center = 'mean', ignore.clusters = NULL, min.cluster.vaf = 0.01, sum.p = 0.05, alpha = 0.05)
  
  y <- transfer.events.to.consensus.trees(y, x, cluster.col.name = 'cluster', event.col.name = 'second')
  
  y <- convert.consensus.tree.clone.to.branch(y, cluster.col = 'cluster', branch.scale = 'sqrt')
  
  
  plot.clonal.models(y, out.dir = 'output', box.plot = TRUE, samples = sample.groups, fancy.boxplot = TRUE, fancy.variant.boxplot.jitter.alpha = 1, fancy.variant.boxplot.jitter.center.color = 'grey50', fancy.variant.boxplot.base_size = 12, fancy.variant.boxplot.plot.margin = 1, fancy.variant.boxplot.vaf.suffix = '.VAF', clone.shape = 'bell', bell.event = FALSE, bell.event.label.color = 'blue', bell.event.label.angle = 60, clone.time.step.scale = 1, bell.curve.step = 2, merged.tree.plot = FALSE, tree.node.label.split.character = NULL, tree.node.shape = 'circle', tree.node.size = 30, tree.node.text.size = 0.5, merged.tree.node.size.scale = 1.25, merged.tree.node.text.size.scale = 2.5, merged.tree.cell.frac.ci = FALSE, merged.tree.clone.as.branch = FALSE, mtcab.event.sep.char = ',', mtcab.branch.text.size = 1, mtcab.branch.width = 0.75, mtcab.node.size = 3, mtcab.node.label.size = 1, mtcab.node.text.size = 1.5, cell.plot = TRUE, num.cells = 100, cell.border.size = 0.25, cell.border.color = 'black', clone.grouping = 'horizontal', scale.monoclonal.cell.frac = TRUE, show.score = FALSE, cell.frac.ci = TRUE, disable.cell.frac = FALSE, overwrite.output = TRUE,  out.format = 'pdf', width = 8, height = 4, panel.widths = c(3, 4, 2))
}






## with VAFs instead of cellular_prevalence
data_for_plot <- PyClone_with_all_mutations %>%
  select(PATIENT_ID, cluster_id, variant_allele_frequency) %>%
  na.omit() %>%
  mutate(variant_allele_frequency = variant_allele_frequency*100)
  
  data_for_plot$cluster_id = data_for_plot$cluster_id+1
  
  # these patients generate unknown errors and cannot be processed with ClonEvol
  data_for_plot = data_for_plot %>%
    filter(PATIENT_ID != "14-00543") %>%
  filter(PATIENT_ID != "14-00667") %>%
  filter(PATIENT_ID != "16-00460") %>%
  filter(PATIENT_ID != "16-00525") %>%
  filter(PATIENT_ID != "TCGA-AB-2907")


dir.create("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Evolution_with_VAFs/")

for (i in 1:length(unique(data_for_plot$PATIENT_ID))) {
  print(i)
  patient <- unique(data_for_plot$PATIENT_ID)[i]
  print(patient)
  dir.create(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Evolution_with_VAFs/", patient))
  setwd(paste0("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Evolution_with_VAFs/", patient))
  height <- length(which(data_for_plot$PATIENT_ID == patient))
  x <- data.frame(cluster = rep(NA, height), bogus.vaf = 0, second.vaf = 0)
  input_data <- filter(data_for_plot, PATIENT_ID == patient)
  
  for (i in 1:height) {
    x$cluster[i] <- input_data$cluster_id[i]
    x[i, 2] <- input_data$variant_allele_frequency[i]
    x[i, 3] <- input_data$variant_allele_frequency[i]
  }
  vaf.col.names <- grep('.vaf', colnames(x), value=T)
  sample.names <- gsub('.vaf', '', vaf.col.names)
  x[, sample.names] <- x[, vaf.col.names]
  vaf.col.names <- sample.names
  sample.groups <- c("bogus", "second")
  names(sample.groups) <- vaf.col.names
  x <- x[order(x$cluster),]
  clone.colors <- NULL
  pdf('box.pdf', width = 3, height = 3, useDingbats = FALSE, title = '')
  plot.variant.clusters(x, cluster.col.name = 'cluster', show.cluster.size = TRUE, cluster.size.text.color = 'blue', vaf.col.names = vaf.col.names, vaf.limits = 100, sample.title.size = 20, box = FALSE, jitter = TRUE, jitter.shape = 1, jitter.color = "black", jitter.center.method = 'median', jitter.alpha = 1)
  dev.off()
  
  x<- x %>%
    group_by(cluster) %>%
    mutate(average_vaf = mean(bogus))
  founding.cluster <- x$cluster[first(which(x$average_vaf == max(x$average_vaf)))]
  
  y <- infer.clonal.models(variants = x, cluster.col.name = 'cluster', vaf.col.names = vaf.col.names, sample.groups = sample.groups, cancer.initiation.model = 'monoclonal', subclonal.test = 'bootstrap',  subclonal.test.model = 'non-parametric', num.boots = 1000, founding.cluster = founding.cluster, cluster.center = 'mean', ignore.clusters = NULL, min.cluster.vaf = 0.01, sum.p = 0.05, alpha = 0.05)
  
  y <- transfer.events.to.consensus.trees(y, x, cluster.col.name = 'cluster', event.col.name = 'second')
  
  y <- convert.consensus.tree.clone.to.branch(y, cluster.col = 'cluster', branch.scale = 'sqrt')
  
  plot.clonal.models(y, out.dir = 'output', box.plot = TRUE, samples = sample.groups, fancy.boxplot = TRUE, fancy.variant.boxplot.jitter.alpha = 1, fancy.variant.boxplot.jitter.center.color = 'grey50', fancy.variant.boxplot.base_size = 12, fancy.variant.boxplot.plot.margin = 1, fancy.variant.boxplot.vaf.suffix = '.VAF', clone.shape = 'bell', bell.event = FALSE, bell.event.label.color = 'blue', bell.event.label.angle = 60, clone.time.step.scale = 1, bell.curve.step = 2, merged.tree.plot = FALSE, tree.node.label.split.character = NULL, tree.node.shape = 'circle', tree.node.size = 30, tree.node.text.size = 0.5, merged.tree.node.size.scale = 1.25, merged.tree.node.text.size.scale = 2.5, merged.tree.cell.frac.ci = FALSE, merged.tree.clone.as.branch = FALSE, mtcab.event.sep.char = ',', mtcab.branch.text.size = 1, mtcab.branch.width = 0.75, mtcab.node.size = 3, mtcab.node.label.size = 1, mtcab.node.text.size = 1.5, cell.plot = TRUE, num.cells = 100, cell.border.size = 0.25, cell.border.color = 'black', clone.grouping = 'horizontal', scale.monoclonal.cell.frac = TRUE, show.score = FALSE, cell.frac.ci = TRUE, disable.cell.frac = FALSE, overwrite.output = TRUE,  out.format = 'pdf', width = 8, height = 4, panel.widths = c(3, 4, 2))
}


## Branched trajectories
# 
# generous assignment (including mixed trajectory models)
# vaf_branched <- c("12-00127", "13-00204", "13-00226", "13-00342", "13-00417", "13-00420", "13-00540", "13-0563", "13-00615", "13-00625", "14-00034", "14-00092", "14-00272", "14-00329", "14-00404", "14-00473", "14-00564", "15-00073", "15-00300", "15-00361", "15-00471", "15-00572", "15-00595", "15-00614", "15-00655", "15-00724", "15-00726", "15-00782", "15-00898", "16-00027", "16-00035", "16-00067", "16-00129", "16-00220", "16-01121", "16-01219", "TCGA-AB-2803", "TCGA-AB-2804", "TCGA-AB-2805", "TCGA-AB-2808", "TCGA-AB-2814", "TCGA-AB-2822", "TCGA-AB-2828", "TCGA-AB-2894", "TCGA-AB-2899", "TCGA-AB-2900", "TCGA-AB-2904", "TCGA-AB-2908", "TCGA-AB-2916", "TCGA-AB-2920", "TCGA-AB-2926", "TCGA-AB-2939", "TCGA-AB-2959", "TCGA-AB-2969", "TCGA-AB-2970", "TCGA-AB-2988", "TCGA-AB-2989", "TCGA-AB-3012")

# # out of curriosity, these are the curated branched samples using cellular prevalence instead of VAF
# cp_branched <- c("12-00127", "14-00272", "14-00404", "15-00073", "15-00595", "15-00724", "16-00035", "16-01151", "TCGA-AB-2802", "TCGA-AB-2808", "TCGA-AB-2817", "TCGA-AB-2830", "TCGA-AB-2897", "TCGA-AB-2924", "TCGA-AB-2965", "TCGA-AB-2969", "TCGA-AB-2978")
# 
# branched <- union(vaf_branched, cp_branched)

# stringent assignment
vaf_branched <- c("12-00127", "13-00540", "13-00563", "13-00625", "14-00272", "14-00329", "14-00404", "14-00473", "14-00564", "14-00632", "15-00073", "15-00361", "15-00572", "15-00595", "15-00724", "15-00726", "15-00898", "16-00035", "16-00220", "16-00525", "16-01151", "16-01219", "TCGA-AB-2803", "TCGA-AB-2808", "TCGA-AB-2814", "TCGA-AB-2828", "TCGA-AB-2897", "TCGA-AB-2900", "TCGA-AB-2916", "TCGA-AB-2959", "TCGA-AB-2969", "TCGA-AB-3012")

branched = vaf_branched

# add a column in data_for_plot for branched/linear evolution
data_for_plot <- PyClone %>%
  select(PATIENT_ID, cluster_id, variant_allele_frequency, overallSurvival, Censor) %>%
  na.omit() %>%
  select(PATIENT_ID, overallSurvival, Censor) %>%
  mutate(Censor = as.numeric(Censor)) %>%
  distinct() %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE))

setwd("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/")

fit <- survfit(Surv(overallSurvival, Censor) ~ is_branched, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, risk.table.col = "strata", ggtheme = theme_bw(), palette = c("black", "blue")) + ggsave("Evolution_survival.pdf", units="in", width=5, height=4, dpi=300) + ggtitle("Evolution Survival")

###########################
# 19. Fisher's Exact Test #
###########################
# 1130 unique gene symbols in TCGA
# 1046 unique gene symbols in BeatAML
# 135 genes overlap
# 1704 unique gene symbols in PyClone

# determine top mutated genes in BeatAML and TCGA
# data_for_plot <- PyClone %>%
#   select(symbol) %>%
#   group_by(symbol) %>%
#   summarize(frequency = length(symbol)) %>%
#   arrange(desc(frequency)) %>%
#   filter(frequency > 5) %>%
#   na.omit()
# genes <- data_for_plot$symbol
# 
# fisher_data <- data.frame(gene_1 = rep(NA, nrow(data_for_plot)^2))
# fisher_data$gene_2 <- NA
# fisher_data$both_present <- 0
# fisher_data$gene_1_present <- 0
# fisher_data$gene_2_present <- 0
# fisher_data$neither_present <- 0
# fisher_data$odds_ratio <- 0
# fisher_data$p_value <- 0
# 
# for (i in 1:length(data_for_plot$symbol)) {
#   gene_1 <- data_for_plot$symbol[i]
#   for (j in 1:length(data_for_plot$symbol)) {
#     gene_2 <- data_for_plot$symbol[j]
#     fisher_data$gene_1[(i-1)*length(data_for_plot$symbol) + j] <- gene_1
#     fisher_data$gene_2[(i-1)*length(data_for_plot$symbol) + j] <- gene_2
#   }
# }
# 
# # exclude repeats
# test <- c()
# fisher_data$include <- NA
# for (i in 1:nrow(fisher_data)) {
#   if (str_c(fisher_data$gene_1[i], fisher_data$gene_2[i], sep = "_") %in% test | str_c(fisher_data$gene_2[i], fisher_data$gene_1[i], sep = "_") %in% test) {
#     fisher_data$include[i] <- FALSE
#   } else {
#     fisher_data$include[i] <- TRUE
#     test <- append(test, str_c(fisher_data$gene_1[i], fisher_data$gene_2[i], sep = "_"))
#     test <- append(test, str_c(fisher_data$gene_2[i], fisher_data$gene_1[i], sep = "_"))
#   }
# }
# fisher_data <- fisher_data[which(fisher_data$include), 1:8]
# 
# data_for_plot <- PyClone %>%
#   select(PATIENT_ID, cluster_id, symbol) %>%
#   na.omit() %>%
#   mutate(unique_id = str_c(PATIENT_ID, cluster_id, sep = "_"))
# 
# for (i in 1:nrow(fisher_data)) {
#   gene_1 <- fisher_data$gene_1[i]
#   gene_2 <- fisher_data$gene_2[i]
#   for (j in 1:length(unique(data_for_plot$unique_id))) {
#     current_id <- unique(data_for_plot$unique_id)[j]
#     if (gene_1 %in% data_for_plot$symbol[which(data_for_plot$unique_id == current_id)]) {
#       if (gene_2 %in% data_for_plot$symbol[which(data_for_plot$unique_id == current_id)]) {
#         fisher_data$both_present[i] <- fisher_data$both_present[i] + 1
#       } else {
#         fisher_data$gene_1_present[i] <- fisher_data$gene_1_present[i] + 1
#       }
#     } else {
#       if (gene_2 %in% data_for_plot$symbol[which(data_for_plot$unique_id == current_id)]) {
#         fisher_data$gene_2_present[i] <- fisher_data$gene_2_present[i] + 1
#       } else {
#         fisher_data$neither_present[i] <- fisher_data$neither_present[i] + 1
#       }
#     }
#   }
# }
# 
# fisher_data$odds_ratio <- (fisher_data$both_present*fisher_data$neither_present)/(fisher_data$gene_1_present*fisher_data$gene_2_present)
# 
# for (i in 1:nrow(fisher_data)) {
#   fishers_test <- matrix(c(fisher_data$both_present[i], fisher_data$gene_2_present[i], fisher_data$gene_1_present[i], fisher_data$neither_present[i]), nrow = 2, dimnames = list(gene_1 = c("present", "absent"), gene_2 = c("present", "absent")))
#   fisher_data$p_value[i] <- fisher.test(fishers_test, alternative = "greater")$p.value
# }
# 
# fisher_data <- fisher_data %>%
#   mutate(keep = ifelse(gene_1 != gene_2, TRUE, FALSE)) %>%
#   filter(keep) %>%
#   mutate(log_odds_ratio = ifelse(odds_ratio == 0, NA, log(odds_ratio, base = 2)))
# 
# data_for_plot <- data.frame(DNMT3A = rep(NA, length(genes)))
# data_for_plot$`FLT3-ITD` = NA
# data_for_plot$NRAS = NA
# data_for_plot$IDH2 = NA
# data_for_plot$TET2 = NA
# data_for_plot$IDH1 = NA
# data_for_plot$`FLT3-TKD` = NA
# data_for_plot$RUNX1 = NA
# data_for_plot$TP53 = NA
# data_for_plot$CEBPA = NA
# data_for_plot$PTPN11 = NA
# data_for_plot$KRAS = NA
# data_for_plot$SRSF2 = NA
# data_for_plot$U2AF1 = NA
# data_for_plot$SMC1A = NA
# data_for_plot$SMC3 = NA
# data_for_plot$WT1 = NA
# data_for_plot$PHF6 = NA
# data_for_plot$RAD21 = NA
# data_for_plot$SF3B1 = NA
# data_for_plot$BCOR = NA
# data_for_plot$GATA2 = NA
# rownames(data_for_plot) <- genes
# for (i in 1:nrow(fisher_data)) {
#   gene_1 = fisher_data$gene_1[i]
#   gene_2 = fisher_data$gene_2[i]
#   data_for_plot[which(rownames(data_for_plot) == gene_1), which(colnames(data_for_plot) == gene_2)] <- fisher_data$log_odds_ratio[i]
# }
# 
# # make NA values the minimum value
# for (i in 1:nrow(data_for_plot)) {
#   for (j in 1:ncol(data_for_plot)) {
#     if (i != j & is.na(data_for_plot[i,j])) {
#       data_for_plot[i,j] <- -0.7049186
#     }
#   }
# }
# 
# p_values <- data.frame(DNMT3A = rep(NA, length(genes)))
# p_values$`FLT3-ITD` = NA
# p_values$NRAS = NA
# p_values$IDH2 = NA
# p_values$TET2 = NA
# p_values$IDH1 = NA
# p_values$`FLT3-TKD` = NA
# p_values$RUNX1 = NA
# p_values$TP53 = NA
# p_values$CEBPA = NA
# p_values$PTPN11 = NA
# p_values$KRAS = NA
# p_values$SRSF2 = NA
# p_values$U2AF1 = NA
# p_values$SMC1A = NA
# p_values$SMC3 = NA
# p_values$WT1 = NA
# p_values$PHF6 = NA
# p_values$RAD21 = NA
# p_values$SF3B1 = NA
# p_values$BCOR = NA
# p_values$GATA2 = NA
# rownames(p_values) <- genes
# for (i in 1:nrow(fisher_data)) {
#   gene_1 = fisher_data$gene_1[i]
#   gene_2 = fisher_data$gene_2[i]
#   p_values[which(rownames(p_values) == gene_1), which(colnames(p_values) == gene_2)] <- fisher_data$p_value[i]
# }
# 
# data_for_plot <- data.matrix(data_for_plot)
# 
# p_values <- data.matrix(p_values)
# corrplot(data_for_plot, method = "circle", type = "upper", is.corr = FALSE, na.label = " ", tl.col = "black", diag = FALSE, addgrid.col = "lightgrey", col = brewer.pal(n = 8, name = "RdBu"), p.mat = p_values, insig = "label_sig", sig.level = c(0.001, 0.01, 0.1), pch.cex = .9, pch.col = "black")

###########
# T tests #
###########

# mutations vs. age
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, total_mutations) %>%
  distinct() %>%
  mutate(age_bin = ifelse(AGE < 60, "< 60", "60 +"))

p_val = round(t.test(data_for_plot$total_mutations[which(data_for_plot$age_bin == "< 60")], data_for_plot$total_mutations[which(data_for_plot$age_bin == "60 +")])$p.val, 3)

data_for_plot %>%
  ggplot(aes(age_bin, total_mutations)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05) + xlab("Age") + ylab("Mutation Burden")  + annotate("text", label = paste("p = ", p_val, sep = ""), x = 1.5, y = 30) + theme_cowplot() + ggsave("Mutation_Burden_vs_Age_Categorical.pdf", units = "in", height = 4, width = 5, dpi = 300)

# clonality vs. age
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, number_of_clones) %>%
  distinct() %>%
  mutate(age_bin = ifelse(AGE < 60, "< 60", "60 +"))

p_val = round( t.test(data_for_plot$number_of_clones[which(data_for_plot$age_bin == "< 60")], data_for_plot$number_of_clones[which(data_for_plot$age_bin == "60 +")])$p.val, 3)

data_for_plot %>%
  ggplot(aes(age_bin, number_of_clones)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05) + xlab("Age") + ylab("Number of Clones")  + annotate("text", label = paste("p = ", p_val, sep = ""), x = 1.5, y = 17) + theme_cowplot() + ggsave("Clonality_vs_Age_Categorical.pdf", units = "in", width = 5, height = 4, dpi = 300)

# age vs. trajectory
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE) %>%
  distinct() %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, "branched", "linear")) %>%
  na.omit()

p_val = round( t.test(data_for_plot$AGE[which(data_for_plot$is_branched == "branched")], data_for_plot$AGE[which(data_for_plot$is_branched == "linear")])$p.val, 3)

data_for_plot %>%
  ggplot(aes(is_branched, AGE)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05) + xlab("Trajectory") + ylab("Age") +  annotate("text", label = paste("p = ", p_val, sep = ""), x = 1.5, y = 90) + theme_cowplot() + ggsave("Age_vs_Trajectory.pdf", units = "in", width = 5, height = 4, dpi = 300)

# Study vs. trajectory
data_for_plot <- PyClone %>%
  select(PATIENT_ID) %>%
  distinct() %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  mutate(study = ifelse(str_detect(PATIENT_ID, "TCGA"), "TCGA", "BeatAML")) %>%
  na.omit()

p_val = round(  t.test(data_for_plot$is_branched[which(data_for_plot$study == "TCGA")], data_for_plot$is_branched[which(data_for_plot$study == "BeatAML")])$p.val, 3)
mean(data_for_plot$is_branched[which(data_for_plot$study == "TCGA")])
# p = 0.3467

########################
# GSEA on Trajectories #
########################
# branched
# data_for_plot <- PyClone %>%
#   select(symbol, PATIENT_ID) %>%
#   mutate(is_branched = ifelse(PATIENT_ID %in% branched, "branched", "linear")) %>%
#   na.omit() %>%
#   filter(is_branched == "branched") %>%
#   group_by(symbol) %>%
#   summarize(frequency = length(symbol)) %>%
#   arrange(desc(frequency)) %>%
#   ungroup() %>%
#   mutate(prevalence = frequency/sum(frequency))
# 
# data_for_plot <- data_for_plot %>%
#   mutate(symbol = ifelse(str_detect(symbol, "FLT3"), "FLT3", symbol))
# 
# 
# data_for_plot$rank <- seq(1, nrow(data_for_plot))
# gene_list <- as.vector(data_for_plot$rank)
# names(gene_list) <- data_for_plot$symbol
# gseResult1 <- EnrichAnalyzer(gene_list, method = "GSEA")
# df1 = gseResult1@result
# df1$logFDR = -log10(df1$pvalue)
# rownames(df1) <- c()
# df1 %>%
#   select(Description, pvalue, logFDR) %>%
#   head(n = 20)
# 
# # linear
# data_for_plot <- PyClone %>%
#   select(symbol, PATIENT_ID) %>%
#   mutate(is_branched = ifelse(PATIENT_ID %in% branched, "branched", "linear")) %>%
#   na.omit() %>%
#   filter(is_branched == "linear") %>%
#   group_by(symbol) %>%
#   summarize(frequency = length(symbol)) %>%
#   arrange(desc(frequency)) %>%
#   ungroup() %>%
#   mutate(prevalence = frequency/sum(frequency))
# 
# data_for_plot <- data_for_plot %>%
#   mutate(symbol = ifelse(str_detect(symbol, "FLT3"), "FLT3", symbol))
# 
# 
# data_for_plot$rank <- seq(1, nrow(data_for_plot))
# gene_list <- as.vector(data_for_plot$rank)
# names(gene_list) <- data_for_plot$symbol
# gseResult1 <- EnrichAnalyzer(gene_list, method = "GSEA")
# df1 = gseResult1@result
# df1$logFDR = -log10(df1$pvalue)
# rownames(df1) <- c()
# df1 %>%
#   select(Description, pvalue, logFDR) %>%
#   head(n = 20)

############################
# 20. Preparing figures ####
############################
if (!require('ggplotify')) install.packages('ggplotify'); library('ggplotify')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('grid')) install.packages('grid'); library('grid')

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

dir.create("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Figure_Graphs/")
setwd("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Figure_Graphs/")
#############################
# Histogram of clone number #
#############################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones) %>%
  distinct() %>%
  filter(!is.na(number_of_clones)) %>%
  mutate(number_of_clones_color = ifelse(number_of_clones <= 5, "1-5", "6+"))

ggplot(data_for_plot, aes(number_of_clones, fill = number_of_clones_color)) + theme_cowplot() + geom_bar(color = "black") + xlab("Number of Clones") + ylab("Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(legend.position = "none") + scale_fill_manual(values = cbPalette) + ggsave("Clone_Histogram.pdf", units = "in", height = 3, width = 4, dpi=300)

############################
# Survival with clone bins #
############################

# PyClone$overallSurvival = PyClone$overallSurvival/365

data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones) %>%
  distinct() %>%
  na.omit() %>%
  mutate(number_of_clones_color = ifelse(number_of_clones <= 5, "1-5", "6+"))


data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
# data_for_plot$two_bins <- NA
# for (i in 1:nrow(data_for_plot)) {
#   if (data_for_plot$number_of_clones[i] <= 3) {
#     data_for_plot$two_bins[i] <- "1-3"
#   }
#   if (data_for_plot$number_of_clones[i] > 3 & data_for_plot$number_of_clones[i] <= 5) {
#     data_for_plot$two_bins[i] <- "3-5"
#   }
#   if (data_for_plot$number_of_clones[i] > 5 & data_for_plot$number_of_clones[i] <= 7) {
#     data_for_plot$two_bins[i] <- "6-7"
#   }
#   if (data_for_plot$number_of_clones[i] > 7) {
#     data_for_plot$two_bins[i] <- "7 +"
#   }
# }

# two bins
fit <- survfit(Surv(overallSurvival, Censor) ~ number_of_clones_color, data = data_for_plot)
ggsurvplot(fit, pval = T, conf.int = FALSE,  xlab = "Years", risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2], cbPalette[3], cbPalette[4]), legend.labs = c("1-5 clones",  "6+ clones"), legend = c(0.8, 0.9), legend.title = "") + ggsave("Survival_Clone_Bins.pdf", units = "in", height = 2.5, width = 3.5, dpi=300)

# # multiple bins
# fit <- survfit(Surv(overallSurvival, Censor) ~ two_bins, data = data_for_plot)
# ggsurvplot(fit, pval = T, conf.int = FALSE,  xlab = "Years", risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2], cbPalette[3], cbPalette[4]), legend.labs = c("1-3 clones", "3-5 clones", "6-7 clones","7+ clones"), legend = c(0.4, 0.9), legend.title = "") + ggsave("Survival_Clone_Bins_granular.pdf", units = "in", height = 3, width = 4, dpi=300)

#############################
# Histogram of mutation bin #
#############################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations) %>%
  distinct() %>%
  filter(!is.na(total_mutations)) %>%
  mutate(number_of_mutations_color = ifelse(total_mutations <= 9, "1-9", "10+"))

ggplot(data_for_plot, aes(total_mutations, fill = number_of_mutations_color)) + theme_cowplot() + geom_bar(color = "black") + xlab("Number of Mutations") + ylab("Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(legend.position = "none") + scale_fill_manual(values = c(cbPalette[3], cbPalette[4])) + ggsave("Mutation_Histogram.pdf", units = "in", height = 3, width = 4, dpi=300)


#############################
# Survival by mutation bins #
#############################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, total_mutations) %>%
  distinct() %>%
  na.omit() %>%
  mutate(number_of_mutations_color = ifelse(total_mutations <= 9, "1-9", "10+"))

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
# data_for_plot$two_bins <- NA
# for (i in 1:nrow(data_for_plot)) {
#   if (data_for_plot$total_mutations[i] <= 5) {
#     data_for_plot$two_bins[i] <- "1-5"
#   }
#   if (data_for_plot$total_mutations[i] > 5 & data_for_plot$total_mutations[i] <= 10) {
#     data_for_plot$two_bins[i] <- "6-10"
#   }
#   if (data_for_plot$total_mutations[i] > 10 & data_for_plot$total_mutations[i] <= 15) {
#     data_for_plot$two_bins[i] <- "10-15"
#   }
#   if (data_for_plot$total_mutations[i] > 15) {
#     data_for_plot$two_bins[i] <- "> 15"
#   }
# }

fit <- survfit(Surv(overallSurvival, Censor) ~ number_of_mutations_color, data = data_for_plot)
ggsurvplot(fit, pval = T, conf.int = FALSE,  xlab = "Years", risk.table = FALSE, palette = c(cbPalette[3], cbPalette[4]),  legend.labs = c("1-9 mutations", "10+ mutations"), legend = c(0.8, 0.9), legend.title = "") + ggsave("Survival_Mutation_Bins.pdf", units = "in", height = 2.5, width = 3.5, dpi=300)


###############################################################################################
# Mutation burden vs. clonality and difference in clonality and mutation burden by tracjetory #
###############################################################################################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, number_of_clones, total_mutations, isDenovo) %>%
  subset(isDenovo == "TRUE") %>% 
  distinct() %>%
  mutate(number_of_clones_color = ifelse(number_of_clones <= 5, "1-5", "6+")) %>%
  mutate(number_of_clones = factor(number_of_clones)) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, "Branched", "Linear")) %>%
  na.omit() 

## with linear and branched stratification
p_kruskal = kruskal.test(total_mutations ~ number_of_clones, data = data_for_plot)$p.val

ggplot(data_for_plot, aes(x = number_of_clones, y = total_mutations)) +
  ylab(label= "# of Driver Mutations") +
  xlab(label = "Number of Clones") +
  theme_cowplot() +
  theme(legend.position=c(0.7, 0.95), legend.title = element_blank())  +
  geom_point(aes(fill = is_branched), subset = .(is_branched == 'Branched'), color = "black", shape = 21, position=position_jitter(0.1), size = 1.5) +
  geom_flat_violin(position = position_nudge(x = 0.3, y = 0),
                   color = "black", fill = "lightgrey",
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
scale_fill_manual(values = c("Linear" = "darkgrey", "Branched" = "#0072B2")) +
  annotate("text", x = 12, y = 1, label = "Kruskal-Wallis, p = 1.6e-33", size = 3) 
  
 ggsave("~/Desktop/MetaAML_results/Figure_5/PyClone/Final_Graphs/Figure_Graphs/Mutations_vs_Clonality_linear_and_branched_raincloud.pdf", units = "in", height = 3, width = 4, dpi=300)







ggplot(data_for_plot, aes(number_of_clones, total_mutations, fill = number_of_clones_color)) + 
  theme_cowplot() + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(data = data_for_plot[which(!data_for_plot$is_branched),], aes(x = number_of_clones, y = total_mutations, alpha = .5), color = "black", width = 0.2) +  
  geom_jitter(data = data_for_plot[which(data_for_plot$is_branched),], aes(x = number_of_clones, y = total_mutations, alpha = .5), color = "red", width = 0.2) + 
  scale_fill_manual(values = cbPalette) +  
  theme(legend.position = "none") + 
  xlab("Number of Clones") + 
  ylab("Number of Mutations") + 
  annotate("text", x = 10, y = 30, label = "Kruskal-Wallis, p = 1.6e-33") + 
  annotate("text", x = 10, y = 28, label = "black = linear evolution") + 
  annotate("text", x = 10, y = 26, label = "red = branched evolution") +
  coord_cartesian(ylim = c(0,30)) + 
  ggsave("Mutations_vs_Clonality_linear_and_branched.pdf", units = "in", height = 3, width = 4, dpi=300)



## without liner and branched stratification
ggplot(data_for_plot, aes(number_of_clones, total_mutations, fill = number_of_clones_color)) + geom_boxplot(outlier.shape = NA) + theme_cowplot() + geom_jitter(data = data_for_plot, color = "black", width = 0.2) + scale_fill_manual(values = cbPalette)  + theme(legend.position = "none") + xlab("Number of Clones") + ylab("Number of Mutations") + annotate("text", x = 10, y = 30, label = "Kruskal-Wallis, p = 1.6e-33") + ggsave("Mutations_vs_Clonality.pdf", units = "in", height = 3, width = 4, dpi=300)

mean(as.numeric(data_for_plot$number_of_clones[which(data_for_plot$is_branched)]), na.rm = TRUE)
# the average number of clones for branched trajectory is 7.62

mean(as.numeric(data_for_plot$number_of_clones[which(!data_for_plot$is_branched)]), na.rm = TRUE)
# the average number of clones for linear trajectory is 4.63

#### t test for difference in number of clones p = 5.284e-08
n_clone_p = signif(t.test(as.numeric(data_for_plot$number_of_clones[which(data_for_plot$is_branched)]), as.numeric(data_for_plot$number_of_clones[which(!data_for_plot$is_branched)]))$p.val, 2)

mean(as.numeric(data_for_plot$total_mutations[which(data_for_plot$is_branched)]), na.rm = TRUE)
# the average number of mutations for branched trajectory is 11.97059

mean(as.numeric(data_for_plot$total_mutations[which(!data_for_plot$is_branched)]), na.rm = TRUE)
# the average number of mutations for linear trajectory is 8.033846

#### t test for difference in number of mutations p = 0.008692
n_mut_p = signif(t.test(as.numeric(data_for_plot$total_mutations[which(data_for_plot$is_branched)]), as.numeric(data_for_plot$total_mutations[which(!data_for_plot$is_branched)]))$p.val, 2)

data_for_plot <- data_for_plot %>%
  mutate(is_branched = ifelse(is_branched, "branched", "linear"))

ggplot(data_for_plot, aes(is_branched, as.numeric(number_of_clones))) + theme_cowplot() + geom_boxplot(aes(fill = is_branched), outlier.shape = NA) + geom_jitter(width = 0.2) + ylab("Number of Clones") + xlab("Evolutionary Trajectory") + annotate("text", x = 1.45, y = 15, label = paste("p = ", n_clone_p)) +  scale_fill_manual(values = c("#0072B2", "darkgrey")) + theme(legend.position = "none") + ggsave("Clonality_by_Trajectory.pdf" , units = "in", height = 3, width = 4, dpi=300)

ggplot(data_for_plot, aes(is_branched, as.numeric(total_mutations))) + theme_cowplot()+ geom_boxplot(aes(fill = is_branched), outlier.shape = NA) + geom_jitter(width = 0.2) + ylab("Number of Mutations") + xlab("Evolutionary Trajectory") + annotate("text", x = 1.45, y = 30, label = paste("p = ", n_mut_p)) + scale_fill_manual(values = c("#0072B2", "darkgrey")) + theme(legend.position = "none") + ggsave("Mutation_Burden_by_Trajectory.pdf" , units = "in", height = 3, width = 4, dpi=300)

#########################################
# Risk by mutation burden and clonality #
#########################################
## By number of mutations
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Risk, total_mutations) %>%
  distinct() %>%
  na.omit()

data_for_plot$Risk[which(data_for_plot$Risk == "IntermediateOrAdverse")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "FavorableOrIntermediate")] <- "Intermediate"
data_for_plot$Risk[which(data_for_plot$Risk == "Good")] <- "Favorable"
data_for_plot$Risk[which(data_for_plot$Risk == "Poor")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "N.D.")] <- "Unknown"
data_for_plot <- data_for_plot %>%
  mutate(Risk = fct_relevel(Risk, "Adverse", "Intermediate", "Favorable", "Unknown"))

ggplot(data_for_plot, aes(total_mutations)) + theme_cowplot() + geom_bar(aes(fill = Risk), position = "stack") + scale_fill_manual(values = cbPalette[5:8]) + xlab("Number of Mutations") + ylab("Number of Patients") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(legend.position = c(.6,.75)) + ggsave("Risk_by_Mutation_Burden.pdf", units = "in", height = 3, width = 4, dpi=300)

## By number of clones
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Risk, number_of_clones) %>%
  distinct() %>%
  na.omit()

data_for_plot$Risk[which(data_for_plot$Risk == "IntermediateOrAdverse")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "FavorableOrIntermediate")] <- "Intermediate"
data_for_plot$Risk[which(data_for_plot$Risk == "Good")] <- "Favorable"
data_for_plot$Risk[which(data_for_plot$Risk == "Poor")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "N.D.")] <- "Unknown"
data_for_plot <- data_for_plot %>%
  mutate(Risk = fct_relevel(Risk, "Adverse", "Intermediate", "Favorable", "Unknown"))

ggplot(data_for_plot, aes(number_of_clones)) + theme_cowplot() + geom_bar(aes(fill = Risk), position = "stack") + scale_fill_manual(values = cbPalette[5:8]) + xlab("Number of Clones") + ylab("Number of Patients") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(legend.position = c(.6,.75)) + ggsave("Risk_by_Clonality.pdf", units = "in", height = 3, width = 4, dpi=300)

############################################################
# Number of clones vs. age and number of mutations vs. age #
############################################################
# Clones vs. age
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, number_of_clones) %>%
  distinct() %>%
  na.omit() %>%
  mutate(number_of_clones_color = ifelse(number_of_clones <= 5, "1-5", "6+")) %>%
  mutate(number_of_clones = factor(number_of_clones))

p_kruskal = signif(kruskal.test(number_of_clones ~ AGE, data = data_for_plot)$p.val, 2)

ggplot(data_for_plot, aes(AGE, number_of_clones, fill = number_of_clones_color)) + theme_cowplot() + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0.2) + xlab("Age") + ylab("Number of Clones") + scale_fill_manual(values = cbPalette) + theme(legend.position = "none") + annotate("text", x = 35, y = 13, label = paste("Kruskal-Wallis, p =",p_kruskal), size=3) + ggsave("Clonality_vs_Age.pdf", units = "in", height = 3, width = 4, dpi=300)


library(ggpubr)
# Grouped Scatter plot with marginal density plots

plot = ggscatterhist(
  data_for_plot, x = "AGE", y = "number_of_clones",
  color = "number_of_clones_color", size = 3,
  palette = cbPalette[1:2],
  margin.params = list(fill = "number_of_clones_color", color = "black", size = 0.2)
)  
p = ggpar(plot,
      xlab = "Age", 
      ylab = "Number of Clones",
      legend.title = "# clones")
p + geom_text(x = 30, y = 15, label = "text") + ggsave("Clonality_vs_Age_Categorical.pdf", units = "in", width = 5, height = 4, dpi = 300)
      

# Mutation burden vs. age
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, total_mutations) %>%
  distinct() %>%
  na.omit() %>%
  mutate(mutation_burden_color = ifelse(total_mutations <= 10, "1-10", "11+")) %>%
  mutate(total_mutations = factor(total_mutations))

p_kruskal = signif(kruskal.test(total_mutations ~ AGE, data = data_for_plot)$p.val, 2)

ggplot(data_for_plot, aes(AGE, total_mutations, fill = mutation_burden_color)) + theme_cowplot() + geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0.2) + xlab("Age") + ylab("Number of Mutations") +  scale_fill_manual(values = cbPalette[3:4]) + theme(legend.position = "none", axis.text.y = element_text(size = 7)) + annotate("text", x = 40, y = 26, label = paste("Kruskal-Wallis, p =",p_kruskal), size=3) + ggsave("Mutation_Burden_vs_Age.pdf", units = "in", height = 3, width = 4, dpi=300)

# Grouped Scatter plot with marginal density plots

plot = ggscatterhist(
  data_for_plot, x = "AGE", y = "total_mutations",
  color = "mutation_burden_color", size = 3,
  palette = cbPalette[3:4],
  margin.params = list(fill = "mutation_burden_color", color = "black", size = 0.2)
)  

p = ggpar(plot,
          xlab = "Age", 
          ylab = "Number of Mutations",
          legend.title = "# mutations")
p + geom_text(x = 30, y = 15, label = "text") + ggsave("Mutation_vs_Age_Categorical.pdf", units = "in", width = 6, height = 5.5, dpi = 300)


############################
# Cox proportional-hazards #
############################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, Censor, overallSurvival, total_mutations, Risk, number_of_clones) %>%
  distinct() %>%
  na.omit() %>%
  mutate(Censor = as.numeric(Censor) + 1) %>%
  mutate(clone_bin = ifelse(number_of_clones <= 5, "1-5", "6+")) %>%
  mutate(mutation_bin = ifelse(total_mutations <= 10, "1-10", "11+")) %>%
  mutate(age_bin = ifelse(AGE < 60, "<60", "60+")) %>%
  mutate(study = ifelse(substr(PATIENT_ID, 1, 1) == "T", "TCGA", "BeatAML")) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, "branched", "linear")) %>%
  mutate(trajectory_and_clone = str_c(clone_bin, is_branched)) %>%
  mutate(trajectory_and_mutation = str_c(mutation_bin, is_branched))

data_for_plot$Risk[which(data_for_plot$Risk == "IntermediateOrAdverse")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "FavorableOrIntermediate")] <- "Intermediate"
data_for_plot$Risk[which(data_for_plot$Risk == "Good")] <- "Favorable"
data_for_plot$Risk[which(data_for_plot$Risk == "Poor")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "N.D.")] <- "Unknown"

coxph(formula = Surv(overallSurvival, Censor) ~ clone_bin, data = data_for_plot)
# 7+ clones hazard ratio = 1.10227, se = 0.14906, p = 0.514

coxph(formula = Surv(overallSurvival, Censor) ~ mutation_bin, data = data_for_plot)
# 11+ mutations hazard ratio = 1.4881, se = 0.1501, p = 0.00811

data_for_plot$is_branched = factor(data_for_plot$is_branched, levels = c("linear", "branched"))

coxph(formula = Surv(overallSurvival, Censor) ~ is_branched, data = data_for_plot)
# branched hazard ratio = 0.4741, se = 0.3616, p = 0.02

coxph(formula = Surv(overallSurvival, Censor) ~ age_bin, data = data_for_plot)
# 60+ hazard ratio = 2.4688, se = 0.1622, p = 2.23e-10

data_for_plot <- data_for_plot %>%
  filter(Risk != "Unknown") %>%
  mutate(Risk = fct_relevel(Risk, "Favorable", "Intermediate", "Adverse"))

coxph(formula = Surv(overallSurvival, Censor) ~ Risk, data = data_for_plot)
# intermediate HR = 2.7835, se = 0.2323, p = 1.05e-05
# advanced HR = 3.5925, se = 0.2446, p = 1.71e-07

coxph(formula = Surv(overallSurvival, Censor) ~ study, data = data_for_plot)
# TCGA hazard ratio = 0.96881, se = 0.15558, = p = 0.839

cox <- data.frame(condition = c("1-6 clones", "7+ clones", "1-10 mutations", "11+ mutations", "Linear Trajectory", "Branched Trajectory", "Age < 60", "Age 60+", "Favorable Risk", "Intermediate Risk", "Adverse Risk", "BeatAML", "TCGA"), value = c(1, 1.10227, 1, 1.4881, 1, 0.4741, 1, 2.4688, 1, 2.7835, 3.5925, 1, 0.96881), se = c(0, 0.14906, 0, 0.1501, 0, 0.3616, 0, 0.1622, 0, 0.2323, 0.2446, 0, 0.15558), color_scheme = c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen"), p = c("", "p = 0.514", "", "p = 0.008", "", "p = 0.02", "", "p = 2.23e-10", "", "p = 1.05e-05", "p = 1.71e-07", "", "p = 0.839")) %>%
  mutate(condition = fct_relevel(condition, "TCGA", "BeatAML", "Adverse Risk", "Intermediate Risk", "Favorable Risk", "Age 60+", "Age < 60", "Branched Trajectory", "Linear Trajectory", "11+ mutations", "1-10 mutations", "7+ clones", "1-6 clones")) %>%
  mutate(color_scheme = fct_relevel(color_scheme, "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen"))

newPalette <- c(cbPalette[1:4], "darkgrey", "#0072B2", "grey", "grey", cbPalette[5:8], "#0073C2FF", "#EFC000FF")

ggplot(cox, aes(condition, value, color = color_scheme, label = p)) + theme_cowplot() + geom_hline(yintercept = c(1, 1), linetype = "dashed") + geom_point(size = 2.5) + coord_flip() + ylim(0, 6) + ylab("Hazard Ratio") + scale_color_manual(values = newPalette) + geom_linerange(aes(ymin = value - se, ymax = value + se)) + theme(legend.position = "none", axis.title.y = element_blank()) + geom_text(size = 2.5, color = "#525252", nudge_y = 1.2, nudge_x = 0) + ggsave("Univariate_Cox.pdf", units = "in", height = 3, width = 4, dpi=300)

## p value table
cox <- cox %>%
  select(p)

library(gridExtra)
pdf("Cox_p.pdf")
grid.table(cox, rows = c())
dev.off()


### Stratifying by trajectory and clonality and trajectory and mutation burden
coxph(formula = Surv(overallSurvival, Censor) ~ trajectory_and_clone, data = data_for_plot)
#                                 coef exp(coef) se(coef)     z     p
# trajectory_and_clone1-5linear  0.7799    2.1813   1.0053 0.776 0.438
# trajectory_and_clone6+branched 0.1272    1.1357   1.0696 0.119 0.905
# trajectory_and_clone6+linear   0.9836    2.6741   1.0074 0.976 0.329

coxph(formula = Surv(overallSurvival, Censor) ~ trajectory_and_mutation, data = data_for_plot)
#                                     coef exp(coef) se(coef)     z       p
# trajectory_and_mutation1-10linear   2.389    10.901    1.005 2.377 0.01744
# trajectory_and_mutation11+branched  2.875    17.725    1.070 2.687 0.00720
# trajectory_and_mutation11+linear    2.692    14.767    1.008 2.672 0.00755


cox <- data.frame(condition = c("Branched 1-6 clones", "Branched 7+ clones", "Linear 1-6 clones", "Linear 7+ clones", "Branched 1-10 mutations", "Branched 11+ mutations", "Linear 1-10 mutations", "Linear 11+ mutations"), value = c(1, 1.1357, 2.1813, 2.6741, 1, 17.725, 10.901, 14.767), se = c(0, 1.0696, 1.0053, 1.0074, 0, 1.070, 1.005, 1.008), color_scheme = c("one", "two", "three", "four", "five", "six", "seven", "eight"), p = c("", "p = 0.905", "p = 0.438", "p = 0.329", "", "p = 0.007", "p = 0.017", "p = 0.008")) %>%
  mutate(condition = fct_relevel(condition, "Linear 11+ mutations", "Linear 1-10 mutations", "Branched 11+ mutations", "Branched 1-10 mutations", "Linear 7+ clones", "Linear 1-6 clones", "Branched 7+ clones", "Branched 1-6 clones")) %>%
  mutate(color_scheme = fct_relevel(color_scheme, "one", "two", "three", "four", "five", "six", "seven", "eight"))

newPalette <- c("#b57d00", "#27546e", cbPalette[1], cbPalette[2], "#005e45", "#aba22e", cbPalette[3], cbPalette[4])

ggplot(cox, aes(condition, value, color = color_scheme, label = p)) + theme_cowplot() + geom_hline(yintercept = c(1, 1), linetype = "dotted") + ylim(0,27) + geom_point(size = 2.5) + coord_flip()  + ylab("Hazard Ratio") + scale_color_manual(values = newPalette) + geom_linerange(aes(ymin = value - se, ymax = value + se)) + theme(legend.position = "none", axis.title.y = element_blank()) + geom_text(size = 2.5, color = "#525252", nudge_y = 5.5, nudge_x = 0) + ggsave("Cox_2.pdf", units = "in", height = 3, width = 4, dpi=300)

## p value table
cox <- cox %>%
  select(p)

library(gridExtra)
pdf("Cox_2_p.pdf")
grid.table(cox, rows = c())
dev.off()

### Cox analysis together
data_for_plot$is_branched = relevel(data_for_plot$is_branched, ref = "linear")

model <- coxph(Surv(overallSurvival, Censor) ~ trajectory_and_mutation + trajectory_and_clone + study + Risk + is_branched + clone_bin + mutation_bin, data = data_for_plot)
ggforest(model) + ggsave("Multivariate_model.pdf", units = "in", height = 5, width = 7, dpi=300)


# just for trajectory modeling
data_for_plot <- PyClone %>%
  select(PATIENT_ID, AGE, Censor, overallSurvival, total_mutations, Risk, number_of_clones) %>%
  distinct() %>%
  na.omit() %>%
  mutate(Censor = as.numeric(Censor) + 1) %>%
  mutate(clone_bin = ifelse(number_of_clones <= 5, "1-5", "6+")) %>%
  mutate(mutation_bin = ifelse(total_mutations <= 10, "1-10", "11+")) %>%
  mutate(age_bin = ifelse(AGE < 60, "<60", "60+")) %>%
  mutate(study = ifelse(substr(PATIENT_ID, 1, 1) == "T", "TCGA", "BeatAML")) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, "branched", "linear")) %>%
  mutate(trajectory_and_clone = str_c(clone_bin, is_branched)) %>%
  mutate(trajectory_and_mutation = str_c(mutation_bin, is_branched))

data_for_plot$Risk[which(data_for_plot$Risk == "IntermediateOrAdverse")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "FavorableOrIntermediate")] <- "Intermediate"
data_for_plot$Risk[which(data_for_plot$Risk == "Good")] <- "Favorable"
data_for_plot$Risk[which(data_for_plot$Risk == "Poor")] <- "Adverse"
data_for_plot$Risk[which(data_for_plot$Risk == "N.D.")] <- "Unknown"

data_for_plot <- data_for_plot %>%
  filter(Risk != "Unknown") %>%
  mutate(Risk = fct_relevel(Risk, "Favorable", "Intermediate", "Adverse"))

model <- coxph(Surv(overallSurvival, Censor) ~ is_branched + age_bin + study + Risk + mutation_bin, data = data_for_plot)
ggforest(model) + ggsave("Multivariate_model.pdf", units = "in", height = 5, width = 7, dpi=300)



##########################
# Survival by trajectory #
##########################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, cluster_id, variant_allele_frequency, overallSurvival, Censor) %>%
  na.omit() %>%
  select(PATIENT_ID, overallSurvival, Censor) %>%
  mutate(Censor = as.numeric(Censor)) %>%
  distinct() %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE))

fit <- survfit(Surv(overallSurvival, Censor) ~ is_branched, data = data_for_plot)
print(fit)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE,  palette = c("darkgrey", "#0072B2"), legend.labs = c("Linear Evolution", "Branched Evolution"), legend = c(0.6, 0.9), legend.title = "") + ggsave("Evolution_survival.pdf", units="in", width=5, height=4, dpi=300) + ggsave("Survival_by_Trajectory.pdf", units = "in", height = 3, width = 4, dpi=300)


################################################
# Survival by trajectory and clone or mutation #
################################################
## Just branched
# clonality
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  filter(is_branched) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= 5) {
    data_for_plot$two_bins[i] <- "1-5"
  } else {
    data_for_plot$two_bins[i] <- "6+"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2]), legend.labs = c("1-5 clones", "6+ clones"), legend = c(0.75, 0.2), legend.title = "") + ggsave("Survival_by_Clonality_Branched_Only.pdf", units = "in", height = 3, width = 4, dpi=300)


# mutation burden
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, total_mutations) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  filter(is_branched) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] <= 10) {
    data_for_plot$two_bins[i] <- "1-10"
  } else {
    data_for_plot$two_bins[i] <- "11+"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, palette = c(cbPalette[3], cbPalette[4]),  legend.labs = c("1-10 mutations", "11+ mutations"),  legend = c(0.75, 0.2), legend.title = "", pval.coord = c(0, 0)) + ggsave("Survival_by_Mutation_Branched_Only.pdf", units = "in", height = 3, width = 4, dpi=300)

### Just linear
# clonality
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  filter(!is_branched) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= 5) {
    data_for_plot$two_bins[i] <- "1-5"
  } else {
    data_for_plot$two_bins[i] <- "6+"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2]), legend.labs = c("1-5 clones", "6+ clones"), legend = c(0.4, 0.9), legend.title = "", pval.coord = c(0, 0)) + ggsave("Survival_by_Clonality_Linear_Only.pdf", units = "in", height = 3, width = 4, dpi=300)

# mutations
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, total_mutations) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  filter(!is_branched) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] <= 10) {
    data_for_plot$two_bins[i] <- "1-10"
  } else {
    data_for_plot$two_bins[i] <- "11+"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, xlab = "Years", risk.table = FALSE, palette = c(cbPalette[3], cbPalette[4]),  legend.labs = c("1-10 mutations", "11+ mutations"), legend = c(0.4, 0.9), legend.title = "", pval.coord = c(0, 0)) + ggsave("Survival_by_Mutation_Linear_Only.pdf", units = "in", height = 3, width = 4, dpi=300)


### stratify by evolution pattern and clones
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= 5) {
    data_for_plot$two_bins[i] <- "1-5"
  } else {
    data_for_plot$two_bins[i] <- "6+"
  }
}

# make a column that concatenates clone bin and evolution pattern
data_for_plot <- data_for_plot %>%
  mutate(clone_ev_bin = str_c(as.character(is_branched), two_bins, sep = "_"))

fit <- survfit(Surv(overallSurvival, Censor) ~ clone_ev_bin, data = data_for_plot)
ggsurvplot(fit, pval = F, conf.int = FALSE, xlab = "Years", risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2], "#b57d00", "#27546e"),  legend = c("right"), legend.title = "", legend.labs = c("Linear 1-5 Clones", "Linear 6+ Clones", "Branched 1-5 Clones", "Branched 6+ Clones")) + ggsave("Survival_by_Clones_and_Trajctory.pdf", units = "in", height = 3, width = 5, dpi=300)


### stratify by evolution pattern and mutation burden
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, total_mutations) %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, TRUE, FALSE)) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] <= 10) {
    data_for_plot$two_bins[i] <- "1-10"
  } else {
    data_for_plot$two_bins[i] <- "11+"
  }
}

# make a column that concatenates evolution bin and evolution pattern
data_for_plot <- data_for_plot %>%
  mutate(mut_ev_bin = str_c(as.character(is_branched), two_bins, sep = "_"))

fit <- survfit(Surv(overallSurvival, Censor) ~ mut_ev_bin, data = data_for_plot)
ggsurvplot(fit, pval = F, conf.int = FALSE, xlab = "Years", risk.table = FALSE, palette = c(cbPalette[3], cbPalette[4], "#005e45", "#aba22e"),  legend = c("right"), legend.title = "", legend.labs = c("Linear 1-10 Mutations", "Linear 11+ Mutations", "Branched 1-10 Mutations", "Branched 11+ Mutations")) + ggsave("Survival_by_Mutations_and_Trajectory.pdf", units = "in", height = 3, width = 6, dpi=300)

########################################################################
# Chi-squared test for goodness-of-fit on median clone number survival #
########################################################################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones) %>%
  distinct() %>%
  na.omit()

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)

fit <- survfit(Surv(overallSurvival, Censor) ~ number_of_clones, data = data_for_plot)
median_survivals <- data.frame("Number_of_Clones" = 1:13, "Median_Survival" = c(705, 517, 462, 325, 499, NA, 359, 931, 800, 365, 782, 140, 600)) %>%
  mutate(Proportion = Median_Survival/sum(Median_Survival, na.rm = TRUE)) %>%
  mutate(Proportion = substr((as.character(Proportion)), 1, 5)) %>%
  rename(`Number of Clones` = Number_of_Clones) %>%
  rename(`Median Survival` = Median_Survival)

library(gridExtra)
pdf("Median_Survivals.pdf")
grid.table(median_survivals, rows = c())
dev.off()


median_survivals <- na.omit(median_survivals)
chisq.test(median_survivals$number_of_clones, median_survivals$median_survival)
## p-value = 0.2329


ggplot(median_survivals, aes(`Number of Clones`, `Median Survival`)) + geom_col(fill = "#999999") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme_half_open()


######################
# Mutation Histogram #
######################
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations) %>%
  distinct() %>%
  filter(!is.na(total_mutations)) %>%
  mutate(total_mutations_color = ifelse(total_mutations <= 10, "1-10", "11+"))

ggplot(data_for_plot, aes(total_mutations, fill = total_mutations_color)) + theme_cowplot() + geom_bar(color = "black") + xlab("Number of Mutations") + ylab("Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  theme(legend.position = "none") + scale_fill_manual(values = cbPalette[3:4]) + ggsave("Mutation_Histogram.pdf", units = "in", height = 3, width = 4, dpi=300)


########################
# Trajectories Barplot #
########################
data_for_plot <- PyClone %>%
  select(PATIENT_ID) %>%
  filter(AGE > 18) %>%
  filter(isDenovo) %>%
  filter(!isRelapse) %>%
  distinct() %>%
  mutate(is_branched = ifelse(PATIENT_ID %in% branched, "branched", "linear"))
# 38 branched and 291 linear

data_for_plot %>%
  ggplot(aes(is_branched)) + theme_cowplot() + geom_bar() + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + xlab("Trajectory") + ylab("Frequency") + annotate("text", label = "n = 38", x = 1, 50) + annotate("text", label = "n = 291", x= 2, y = 300) + ggsave("Trajectories_barplot.pdf", units = "in", height = 3, width = 4, dpi = 305)

#######################
# Mutations per Clone #
#######################
data_for_plot <- PyClone %>%
  filter(AGE > 18) %>%
  filter(isDenovo) %>%
  filter(!isRelapse) %>%
  select(PATIENT_ID, Censor, overallSurvival, total_mutations, number_of_clones) %>%
  distinct() %>%
  mutate(mutations_per_clone = total_mutations/number_of_clones) %>%
  na.omit()

# histogram
data_for_plot %>%
  ggplot(aes(mutations_per_clone)) + geom_histogram(color = "black", fill = "lightblue", binwidth = 1)+ xlab("Mutations Per Clone") + ylab("Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme_cowplot() + ggsave("Mutations_per_Clone_Histogram.pdf", units = "in", width = 4, height = 3, dpi = 305)

# stratify by mean
mean_mutations_per_clone <- mean(data_for_plot$mutations_per_clone)

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_mutation_per_clone_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$mutations_per_clone[i] <= mean_mutations_per_clone) {
    data_for_plot$two_mutation_per_clone_bins[i] <- "Low"
  } else {
    data_for_plot$two_mutation_per_clone_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_mutation_per_clone_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c("black", "blue"),  legend = c(0.4, 0.9), legend.title = "", legend.labs = c("High Mutation per Clone", "Low Mutation per Clone")) + ggsave("Survival_by_Mutation_Per_Clone_Mean.pdf", units = "in", height = 4.5, width = 5, dpi=300)

# stratify by median
median_mutations_per_clone <- median(data_for_plot$mutations_per_clone)

data_for_plot$two_mutation_per_clone_bins_median <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$mutations_per_clone[i] <= median_mutations_per_clone) {
    data_for_plot$two_mutation_per_clone_bins_median[i] <- "Low"
  } else {
    data_for_plot$two_mutation_per_clone_bins_median[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_mutation_per_clone_bins_median, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c("darkred", "darkgrey"), legend = c(0.5, 0.9), legend.title = "", legend.labs = c("High Mutation\nburden per Clone\n", "Low Mutation\nburden per Clone")) + ggsave("Survival_by_Mutation_Per_Clone_Median.pdf", units = "in", height = 3, width = 4, dpi=300)

###############################
# Maximum Mutations per Clone #
###############################
data_for_plot <- PyClone %>%
  filter(AGE > 18) %>%
  filter(isDenovo) %>%
  filter(!isRelapse) %>%
  select(PATIENT_ID, cluster_id, size, Censor, overallSurvival) %>%
  na.omit() %>%
  group_by(PATIENT_ID) %>%
  filter(size == max(size)) %>%
  select(PATIENT_ID, size, Censor, overallSurvival) %>%
  distinct()

# histogram
data_for_plot %>%
  ggplot(aes(size)) + geom_histogram(color = "black", fill = "lightblue", binwidth = 1) + xlab("Maximum Mutations Per Clone") + ylab("Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme_cowplot() + ggsave("Max_Mutations_per_Clone_Histogram.pdf", units = "in", width = 4, height = 3, dpi = 305)

# stratify by mean
mean_max_mutations_per_clone <- mean(data_for_plot$size)

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$two_max_mutation_per_clone_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$size[i] <= mean_max_mutations_per_clone) {
    data_for_plot$two_max_mutation_per_clone_bins[i] <- "Low"
  } else {
    data_for_plot$two_max_mutation_per_clone_bins[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_max_mutation_per_clone_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c("black", "blue"),  legend = c(0.4, 0.9), legend.title = "", legend.labs = c("High Max Mutation per Clone", "Low Max Mutation per Clone")) + ggsave("Survival_by_Max_Mutation_Per_Clone_Mean.pdf", units = "in", height = 4.5, width = 5, dpi=300)

# stratify by median
median_max_mutations_per_clone <- median(data_for_plot$size)

data_for_plot$two_max_mutation_per_clone_bins_median <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$size[i] <= median_max_mutations_per_clone) {
    data_for_plot$two_max_mutation_per_clone_bins_median[i] <- "Low"
  } else {
    data_for_plot$two_max_mutation_per_clone_bins_median[i] <- "High"
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ two_max_mutation_per_clone_bins_median, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c("black", "blue"), ggtheme = theme_half_open(), legend = c(0.4, 0.9), legend.title = "", legend.labs = c("High Max Mutation per Clone", "Low Max Mutation per Clone")) + ggsave("Survival_by_Max_Mutation_Per_Clone_Median.pdf", units = "in", height = 4.5, width = 5, dpi=300)

#################
# 21. Quartiles #
#################
dir.create("~/Desktop/Figure_5/PyClone/Final_Graphs/Quartiles")
setwd("~/Desktop/Figure_5/PyClone/Final_Graphs/Quartiles")

## Number of clones
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, number_of_clones) %>%
  distinct() %>%
  filter(!is.na(Censor)) %>%
  filter(!is.na(overallSurvival)) %>%
  filter(!is.na(number_of_clones))

first_quartile <- quantile(data_for_plot$number_of_clones, 0.25)
third_quartile <- quantile(data_for_plot$number_of_clones, 0.75)

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$clone_quartile_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$number_of_clones[i] <= first_quartile) {
    data_for_plot$clone_quartile_bins[i] <- "low"
  } else if (data_for_plot$number_of_clones[i] >= third_quartile) {
    data_for_plot$clone_quartile_bins[i] <- "high"
  } else {
    data_for_plot$clone_quartile_bins[i] <- NA
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ clone_quartile_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2]), legend = c(0.4, 0.9), legend.labs = c("Upper Quartile Clonality", "Lower Quartile Clonality"), legend.title = "") + ggsave("Clonality_quartiles.pdf", units = "in", height = 3, width = 4, dpi=300)

## Mutation Burden
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, total_mutations) %>%
  distinct() %>%
  filter(!is.na(Censor)) %>%
  filter(!is.na(overallSurvival)) %>%
  filter(!is.na(total_mutations))

first_quartile <- quantile(data_for_plot$total_mutations, 0.25)
third_quartile <- quantile(data_for_plot$total_mutations, 0.75)

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$mutation_quartile_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$total_mutations[i] <= first_quartile) {
    data_for_plot$mutation_quartile_bins[i] <- "low"
  } else if (data_for_plot$total_mutations[i] >= third_quartile) {
    data_for_plot$mutation_quartile_bins[i] <- "high"
  } else {
    data_for_plot$mutation_quartile_bins[i] <- NA
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ mutation_quartile_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c(cbPalette[1], cbPalette[2]), legend = c(0.2, 0.9), legend.labs = c("Upper Quartile Mutation Burden", "Lower Quartile Mutation Burden"), legend.title = "") + ggsave("Mutation_Burden_Quartiles.pdf", units = "in", height = 3, width = 4, dpi=300)


## SDI
data_for_plot <- PyClone %>%
  select(PATIENT_ID, Censor, overallSurvival, SDI) %>%
  distinct() %>%
  filter(!is.na(Censor)) %>%
  filter(!is.na(overallSurvival)) %>%
  filter(!is.na(SDI))

first_quartile <- quantile(data_for_plot$SDI, 0.25)
third_quartile <- quantile(data_for_plot$SDI, 0.75)

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$SDI_quartile_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$SDI[i] <= first_quartile) {
    data_for_plot$SDI_quartile_bins[i] <- "low"
  } else if (data_for_plot$SDI[i] >= third_quartile) {
    data_for_plot$SDI_quartile_bins[i] <- "high"
  } else {
    data_for_plot$SDI_quartile_bins[i] <- NA
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ SDI_quartile_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c(cbPalette[5], cbPalette[6]),  legend = c(0.4, 0.9), legend.labs = c("Upper Quartile SDI", "Lower Quartile SDI"), legend.title = "") + ggsave("SDI_Quartiles.pdf", units = "in", height = 3, width = 4, dpi=300)

## MATH
data_for_plot <- PyClone %>%
  select(PATIENT_ID, total_mutations, number_of_clones, variant_allele_frequency, overallSurvival, Censor, SDI, AGE) %>%
  filter(total_mutations > 1) %>%
  filter(!is.na(Censor)) %>%
  filter(!is.na(overallSurvival)) %>%
  group_by(PATIENT_ID) %>%
  mutate(MATH_Score = 100*mad(variant_allele_frequency)/median(variant_allele_frequency)) %>%
  select(PATIENT_ID, total_mutations, number_of_clones, overallSurvival, Censor, MATH_Score, SDI, AGE) %>%
  filter(complete.cases(MATH_Score)) %>%
  distinct()

first_quartile <- quantile(data_for_plot$MATH_Score, 0.25)
third_quartile <- quantile(data_for_plot$MATH_Score, 0.75)

data_for_plot$overallSurvival <- as.numeric(data_for_plot$overallSurvival)
data_for_plot$Censor <- as.integer(data_for_plot$Censor)
data_for_plot$MATH_quartile_bins <- NA
for (i in 1:nrow(data_for_plot)) {
  if (data_for_plot$MATH_Score[i] <= first_quartile) {
    data_for_plot$MATH_quartile_bins[i] <- "low"
  } else if (data_for_plot$MATH_Score[i] >= third_quartile) {
    data_for_plot$MATH_quartile_bins[i] <- "high"
  } else {
    data_for_plot$MATH_quartile_bins[i] <- NA
  }
}

fit <- survfit(Surv(overallSurvival, Censor) ~ MATH_quartile_bins, data = data_for_plot)
ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE, palette = c(cbPalette[5], cbPalette[6]), legend = c(0.4, 0.9), legend.labs = c("Upper Quartile MATH", "Lower Quartile MATH"), legend.title = "") + ggsave("MATH_Quartiles.pdf", units = "in", height = 3, width = 4, dpi=300)


