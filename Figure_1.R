# ========================================================================================================================================== #
# Figure_1.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 09/01/2020
# Description: Curation of multiple sequencing studies in AML into one dataframe for meta analysis of mutational and clinical correlations. This script will pull published data from AML sequencing studies, homogenize annotations and coding where possible, and plot an UpSet plot and oncoprint of the cohort as seen in Figure 1 of the manuscript Benard et al. "Clonal architecture and variant allele frequency correlate with clinical outcomes and drug response in acute myeloid leukemia".
# ========================================================================================================================================== #

# =================== #
# Load libraries ####
# =================== #
if (!require('scales')) install.packages('scales'); library('scales')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('UpSetR')) install.packages('UpSetR'); library('UpSetR')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('RCurl')) install.packages('RCurl'); library('RCurl')
if (!require('reshape')) install.packages('reshape'); library('reshape')

# =============================== #
# make directory dynamic per user
# =============================== #
dir.create("~/Desktop/MetaAML_results/raw_data")
setwd("~/Desktop/MetaAML_results/raw_data")

# ========================================================================================== #
# Using the different datasets, create a single matrix containing all patients and variables
# ========================================================================================== #

# ========= #
# TCGA ####
# ========= #

# download all AML TCGA data
download.file("http://download.cbioportal.org/laml_tcga_pub.tar.gz", destfile = "~/Desktop/MetaAML_results/raw_data/laml_tcga_pub.tar.gz")
untar("~/Desktop/MetaAML_results/raw_data/laml_tcga_pub.tar.gz", files = c("data_mutations_extended.txt", "data_clinical_patient.txt", "data_CNA.txt"), exdir = "raw_data/")

# additional mutation data
download.file("https://api.gdc.cancer.gov/data/0d8851d7-1af0-4054-a527-5db763138400", destfile = "~/Desktop/MetaAML_results/raw_data/supplementalTable06.tsv")

# read in TCGA variants data
TCGAL_variants_raw <-read.table("~/Desktop/MetaAML_results/raw_data/SupplementalTable06.tsv", header = T, sep = "\t", stringsAsFactors = F, fill = T)

# subset to useful columns
variants <- TCGAL_variants_raw %>%
  dplyr::select("UPN", "TCGA_id", "gene_name", "type", "trv_type", "TumorVAF", "amino_acid_change")

# remove silent mutations and rows containing "-"
variants <- variants[!grepl("silent", variants$trv_type),]
variants <- variants[!grepl("intronic", variants$trv_type),]
variants <- variants[!grepl("-", variants$gene_name),]

# read in TCGA variants for analysis data. This list of mutations is the high quality calls that they use in their analysis and papers
TCGA_variants <-read.table("~/Desktop/MetaAML_results/raw_data/data_mutations_extended.txt", header = T, sep = "\t", stringsAsFactors = F)

# subset to mutations used for analysis from TCGA
variants <- setDT(variants)[gene_name %chin% TCGA_variants$Hugo_Symbol]

# subset to useful columns
variants_2 <- variants %>%
  dplyr::select("gene_name")

# count the number of mutations per gene
mut_table <- aggregate(data.frame(count = variants_2), list(value = variants_2$gene_name), length)
mut_table <- dplyr::select(mut_table, "value", "gene_name")
mut_table = mut_table[-1,]
colnames(mut_table) <- c("gene_name", "gene_mut_count")

variants <- dplyr::left_join(variants, mut_table, by = "gene_name")

# count the number of mutations per patient and add to recurrent genes
mut_table_pts <- aggregate(data.frame(count = variants), list(value = variants$UPN), length)
mut_table_pts <- dplyr::select(mut_table_pts, "value", "count.UPN")
colnames(mut_table_pts) <- c("UPN", "sample_mut_count")

# combine for final data frame
mut_table_final_TCGA <- dplyr::left_join(variants, mut_table_pts, by = "UPN")

# change VAF to match other datasets
# mut_table_final_TCGA$TumorVAF <- (mut_table_final_TCGA$TumorVAF/100)

# make column identifying subset type (i.e. de novo) and cohort
mut_table_final_TCGA$Subset <- "de_novo"
mut_table_final_TCGA$Cohort <- "TCGA"


# ============ #
# Beat AML ####
# ============ #
# Beat AML (Tyner et al.)
# cohort and clinical data
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx")

# read in sample types
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# subset to patients with exome sequencing
pt_subset_2 <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$exomeSeq == "y")

# create a new column to annotate the type of sample and then loop through the clinical annotations and decide if the sample is de novo, secondary, or relapse. This is necessary because some samples are relapse samples from an intially de novo, therapy, or secondary patient
pt_subset_2$Subset <- NA

for(i in 1:nrow(pt_subset_2)){
  if (pt_subset_2$isDenovo[i] == "TRUE" & pt_subset_2$isRelapse[i] == "FALSE") {
    pt_subset_2$Subset[i] <- "de_novo"
  }
  if (pt_subset_2$isDenovo[i] == "TRUE" & pt_subset_2$isRelapse[i] == "TRUE") {
    pt_subset_2$Subset[i] <- "relapse"
  }
  if (pt_subset_2$isRelapse[i] == "TRUE" & pt_subset_2$isTransformed[i] == "FALSE") {
    pt_subset_2$Subset[i] <- "relapse"
  }
  if (pt_subset_2$isRelapse[i] == "TRUE" & pt_subset_2$isTransformed[i] == "TRUE") {
    pt_subset_2$Subset[i] <- "relapse"
  }
  if (pt_subset_2$isTransformed[i] == "TRUE" & pt_subset_2$isRelapse[i] == "FALSE") {
    pt_subset_2$Subset[i] <- "transformed"
  }
  if (pt_subset_2$isDenovo[i] == "FALSE" & pt_subset_2$isRelapse[i] == "FALSE" & pt_subset_2$isTransformed[i] == "FALSE") {
    pt_subset_2$Subset[i] <- "other"
  }
}

pt_subset_2 <- dplyr::select(pt_subset_2, LabId, Subset, PatientId, specimenType)
colnames(pt_subset_2)[1] <- "labId"


# Mutations
# read in BeatAML variants for analysis
BeatAML_variants <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7)

# extract useful columns
BeatAML_variants_sub <- BeatAML_variants %>%
  dplyr::dplyr::select("labId", "symbol", "t_vaf", "variant_class", "short_aa_change")

# make sure column names are consistent across datasets
colnames(BeatAML_variants_sub)[4:5] <- c("variant_type", "amino_acid_change")

# remove duplicate calls
BeatAML_variants_sub <- BeatAML_variants_sub %>%
  distinct(labId, symbol, t_vaf, .keep_all = T)

# filter to variants present in the patient group and ensure no duplicate pts are represented
BeatAML_variants_sub <- setDT(BeatAML_variants_sub)[labId %chin% pt_subset_2$labId]
BeatAML_variants_sub <- unique(BeatAML_variants_sub)
colnames(BeatAML_variants_sub)[3] <- "VAF"
# combined <- BeatAML_variants_sub
BeatAML_variants_sub$VAF <- as.numeric(BeatAML_variants_sub$VAF)
BeatAML_variants_sub$VAF <- BeatAML_variants_sub$VAF*100


colnames(BeatAML_variants_sub)[4] <- c("variant_type")

combined <- left_join(BeatAML_variants_sub, pt_subset_2, by = "labId")
combined$Cohort <- "Tyner"


# combine the BeatAML and TCGA mutation files together to find mutation frequencies and co-occurences
# dplyr::select useful columns
mut_table_final_TCGA_sub <- mut_table_final_TCGA[,c(2,3,6,4,7,10,1,11)]

# homogenize variant type annotations between the two datasets
for (i in 1:nrow(mut_table_final_TCGA_sub)) {
  if(mut_table_final_TCGA_sub$type[i] == "DEL"){
    mut_table_final_TCGA_sub$type[i] <- "Deletion"
  }
  if(mut_table_final_TCGA_sub$type[i] == "INS"){
    mut_table_final_TCGA_sub$type[i] <- "Insertion"
  }
  if(mut_table_final_TCGA_sub$type[i] == "SNP"){
    mut_table_final_TCGA_sub$type[i] <- "SNV"
  }
}

# make sure the column headers are the same for rowbind
cnames <- colnames(combined)
mut_table_final_TCGA_sub$specimenType <- "bone_marrow"
mut_table_final_TCGA_sub <- mut_table_final_TCGA_sub[c(1:7,9,8)]
colnames(mut_table_final_TCGA_sub) <- cnames

# remove duplicate rows in both column prior to merging
mut_table_final_TCGA_sub <- unique(mut_table_final_TCGA_sub)
colnames(mut_table_final_TCGA_sub)[3] <- "VAF"
mut_table_final_TCGA_sub$VAF <- as.numeric(mut_table_final_TCGA_sub$VAF)

# ======================= #
# Combine dataframes ####
# ======================= #
# bind the TCGA and BeatAML per patient mutation dataframes together
combined1 <- rbind.fill(combined, mut_table_final_TCGA_sub)


# ============ #
# Majeti  ####
# ============ #
# ===================================================================================================================================== #
# This bit of code loops through individual txt files for mutation calls for each patient and combines them together into one dataframe.
# The user will need to have downloaded the 'Archive' folder of results from the Majeti lab and have it in the appropriate file path.
# ===================================================================================================================================== $
# load list of file names and rename column header
files_list <- as.data.frame(list.files(path="~/Desktop/MetaAML_results/raw_data/Archive/", pattern="*.txt", full.names=F, recursive=FALSE))

colnames(files_list)[1] <- "file_path"

# create empty list to populate results into
results_list <- list()
z <- 1

# loop through the data frame of individual file names and extract the individual mutations and VAFs
for(i in 1:nrow(files_list)){
  # define the path to the individual file for each itteration
  file_name <- as.character(files_list$file_path[i])
  results <- read.table(paste("~/Desktop/MetaAML_results/raw_data/Archive/",file_name, sep = ""), header=F, fill = T, sep = "\t") # load file

  # strip the patient ID from the file hame
  pt_id <- gsub("\\_.*","", file_name) %>%
    as.character(basename(file_name))

  # use the unique nature of these results files to strip VAFS using the % symbol amd INDEL frequencies using the "." symbol
  vaf_hits <- results[grep("%", results$V10), ]
  indel_hits <- results[grepl("\\.", results$V12), ]

  n_hits_vaf <- as.numeric(nrow(vaf_hits))

  # create new data frame to populate VAF results
  temp_dat1 <- data.frame(matrix(NA, nrow = n_hits_vaf, ncol = 5))
  names(temp_dat1) <- c("Sample", "Mutation", "VAF", "SNV_or_Indel", "amino_acid_change")

  temp_dat1[,1] <- pt_id
  temp_dat1[,2] <- vaf_hits$V8
  temp_dat1[,3] <- vaf_hits$V10
  temp_dat1[,4] <- "SNV"
  temp_dat1[,5] <- vaf_hits$V9

  # create new data frame to populate INDEL results
  n_hits_indel <- as.numeric(nrow(indel_hits))

  if(n_hits_indel > 0){
    temp_dat2 <- data.frame(matrix(NA, nrow = n_hits_indel, ncol = 5))
    names(temp_dat2) <- c("Sample", "Mutation", "VAF", "SNV_or_Indel", "amino_acid_change")

    temp_dat2[,1] <- pt_id
    temp_dat2[,2] <- indel_hits$V4
    temp_dat2[,3] <- indel_hits$V12
    temp_dat2[,4] <- "INDEL"
    temp_dat2[,5] <- indel_hits$V5

    # bind the datasets together
    temp_dat_final <- rbind(temp_dat1, temp_dat2)

    results_list[[z]] <- temp_dat_final
    z <- z + 1
  } else {
    results_list[[z]] <- temp_dat1

    z <- z + 1
  }
}

# bind the list of lists together to generate the final data frame
temp_final = as.data.frame(do.call(rbind, results_list))

# read in the clinical data for the Stanford samples
clin_annotations <- read_excel("~/Desktop/MetaAML_results/raw_data/Archive/160104_RM_AML_Sample_Info.xlsx")

# dplyr::select columns of interest
labels <- dplyr::select(clin_annotations, Sample, `Disease Status`, `% Blasts`)
colnames(labels)[1] <- "labId"

# homogenize the sample type annotations
for(i in 1:nrow(labels)){
  if(labels$`Disease Status`[i] == "De Novo" | labels$`Disease Status`[i] == "De novo"){
    labels$`Disease Status`[i] <- "de_novo"
  }
  if(labels$`Disease Status`[i] == "Relapsed" | labels$`Disease Status`[i] == "Refractory"){
    labels$`Disease Status`[i] <- "relapse"
  }
  if(labels$`Disease Status`[i] == "Secondary"){
    labels$`Disease Status`[i] <- "secondary"
  }
  if(labels$`Disease Status`[i] == "?"){
    labels$`Disease Status`[i] <- "other"
  }
  if(labels$`Disease Status`[i] == "Unknown"){
    labels$`Disease Status`[i] <- "other"
  }
}

colnames(labels)[2] <- "Subset"

# pull the FLT3-ITD vs TKD calls for each patient
clin_annotations <- dplyr::select(clin_annotations, Sample, `Clonal Mutations (genotyped by Ryan)`)
flt3_itd <- clin_annotations[grepl("FLT3-ITD", clin_annotations$`Clonal Mutations (genotyped by Ryan)`), ]
colnames(flt3_itd)[2] <- "Mutation"
flt3_itd$Mutation <- "FLT3"
flt3_itd$variant_type <- "ITD"

# add a new row for the double FLT3 mutation annotation
n = 1
flt3_list = list()

for(i in 1:nrow(flt3_itd)){
  pt = subset(flt3_itd, Sample = flt3_itd$Sample[i])
  pt = pt[1,]
  pt$variant_type = "ITD"
  flt3_list[[z]] <- pt
  z <- z + 1
}
final_flt3_list = as.data.frame(do.call(rbind, flt3_list))

temp_final_flt3 <- rbind.fill(temp_final, final_flt3_list)

# label the mutation as FLT3-ITD or FLT3-TKD
for(i in 1:nrow(temp_final_flt3)){
  if(!is.na(temp_final_flt3$variant_type[i])){
    if(temp_final_flt3$variant_type[i] == "ITD"){
      temp_final_flt3$SNV_or_Indel[i] <- "ITD"
    }
  }
}

# # remove unnecessary column
temp_final_flt3$variant_type <- NULL

colnames(temp_final_flt3)[1:5] <- c("labId", "symbol", "VAF", "variant_type", "amino_acid_change")

# ========================================================================================================================================================= #
# seems like there are rare mutations in the stanford samples that are not recurrent in the other studies. For now, filter to only mutations found in TCGA
# ========================================================================================================================================================= #
# temp_final_flt3 <- setDT(temp_final_flt3)[symbol %chin% variants$gene_name]

# remove the % symbol and normalize the VAF calls
temp_final_flt3$VAF <- sub("%$", "", temp_final_flt3$VAF)
temp_final_flt3$VAF <- as.numeric(as.character(temp_final_flt3$VAF))
# temp_final_flt3$VAF <- (temp_final_flt3$VAF/100)

# add the sample type column
temp_final_flt3 <- left_join(temp_final_flt3, labels, by = "labId")
temp_final_flt3$Cohort <- "Stanford"

# ======================= #
# Combine dataframes ####
# ======================= #
temp_final_flt3$PatientId = temp_final_flt3$labId

combined2 <- rbind.fill(combined1, temp_final_flt3)

combined2$`% Blasts` = NULL


# ======================== #
# Papaemmanuil et al. ####
# ======================== #
# clinical annotations
download.file("https://github.com/gerstung-lab/AML-multistage/blob/master/data/AMLSG_Clinical_Anon.RData?raw=true", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")

load("~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")

# Papaemmanuil et al. mutations
x <- getURL("https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/AMLSG_Genetic.txt")

kb_data_muts <- read.table(text = x, header = T, stringsAsFactors = FALSE)

# use flt3 calls to annotation mutations according to type
kb_data_muts_sub <- dplyr::select(kb_data_muts, SAMPLE_NAME, VAF, VARIANT_TYPE, GENE, AA_CHANGE, Study_Variant_ID)

clinicalData2 <- read.table("~/Desktop/MetaAML_results/raw_data/AML_knowledge_bank_data_clinical.txt", stringsAsFactors = F)

# save data for survival analysis
Multistage_survival <- clinicalData2

clin_annotations2 <- dplyr::select(clinicalData2, PDID, TypeAML, FLT3_ITD, FLT3_TKD)
colnames(clin_annotations2)[1] <- "SAMPLE_NAME"

kb_data_flt3 <- left_join(kb_data_muts_sub, clin_annotations2,  by = "SAMPLE_NAME")

kb_data_flt3[is.na(kb_data_flt3)] <- "empty"


## for now, don't specify the number of ITDs reported in the sample. If there are both TKD and ITD in the same sample, default to annotated VARIANT_TYPE"
# rename FLT3 variant type based on the mutation type
for(i in 1:nrow(kb_data_flt3)){
  if(kb_data_flt3$FLT3_TKD[i] == 1 & kb_data_flt3$FLT3_ITD[i] == 0){
    kb_data_flt3$VARIANT_TYPE[i] <- "SNV"
  }
  if(kb_data_flt3$FLT3_ITD[i] == 1 & kb_data_flt3$GENE[i] == "FLT3"){
    kb_data_flt3$VARIANT_TYPE[i] <- "ITD"
  }
  if(kb_data_flt3$VARIANT_TYPE[i] == "Sub"){
    kb_data_flt3$VARIANT_TYPE[i] <- "SNV"
  }
  if(kb_data_flt3$TypeAML[i] == "empty"){
    kb_data_flt3$TypeAML[i] <- "oAML"
  }
}

names <- colnames(kb_data_flt3)

joined_muts_clinical <- left_join(kb_data_muts_sub, kb_data_flt3,  by = "SAMPLE_NAME")

joined_muts_clinical <- joined_muts_clinical[,c(1:6,12:14)]
colnames(joined_muts_clinical) <- c(names)


# rename the varient type for homogeneity with the other data sets
for(i in 1:nrow(joined_muts_clinical)){

  if(joined_muts_clinical$VARIANT_TYPE[i] == "D" & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "Deletion"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "I" & joined_muts_clinical$FLT3_ITD[i] == 1 & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "ID" & joined_muts_clinical$FLT3_ITD[i] == 1 & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "tandem_duplication" & joined_muts_clinical$FLT3_ITD[i] == 1 & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "D" & joined_muts_clinical$GENE[i] != "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "Deletion"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "I" & joined_muts_clinical$GENE[i] != "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "Insertion"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "ID" & joined_muts_clinical$GENE[i] != "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "INDEL"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "Sub"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "SNV"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "PTD"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
}

joined_muts_clinical[,8:9] <- NULL
joined_muts_clinical <- unique(joined_muts_clinical)

joined_muts_clinical <- dplyr::select(joined_muts_clinical, SAMPLE_NAME, GENE, VARIANT_TYPE, VAF, TypeAML, AA_CHANGE, Study_Variant_ID)
colnames(joined_muts_clinical)[1:7] <- c("labId", "symbol", "variant_type", "VAF", "Subset", "amino_acid_change", "Study_Variant_ID")

for(i in 1:nrow(joined_muts_clinical)){
  if(joined_muts_clinical$Subset[i] == "AML"){
    joined_muts_clinical$Subset[i] <- "de_novo"
  }
  if(joined_muts_clinical$Subset[i] == "sAML"){
    joined_muts_clinical$Subset[i] <- "transformed"
  }
  if(joined_muts_clinical$Subset[i] == "rAML"){
    joined_muts_clinical$Subset[i] <- "relapse"
  }
  if(joined_muts_clinical$Subset[i] == "tAML"){
    joined_muts_clinical$Subset[i] <- "therapy"
  }
  if(joined_muts_clinical$Subset[i] == "oAML"){
    joined_muts_clinical$Subset[i] <- "other"
  }
}

# normalize the VAF to decimal format to match other datasets
joined_muts_clinical$VAF <- as.numeric(joined_muts_clinical$VAF, na.rm = T)
joined_muts_clinical$VAF <- (joined_muts_clinical$VAF/100)

# add cohort annotation
joined_muts_clinical$Cohort <- "Papaemmanuil"

###
##
# now, modifying some supplimental code from Gerstung et al. 2016, perform copy number correction on VAFs using avaliable karyotype data and append a column to the mutation table

# 1.3.2.1 Clinical data
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
                Cytogenetics = clinicalData[,grep("^(t_)|(inv)|(abn)|(plus)|(minus)|(mono)|(complex)",colnames(clinicalData))])
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

# download the karyotype data
download.file("https://github.com/gerstung-lab/AML-multistage/raw/master/data/AMLSG_Karyotypes.txt", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Karyotypes.txt")

# here is the code that accounts for broad copy number status correction of VAFs based on karyotype
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
mutationData$CN_status =mutationData$CN_status %>% replace_na(2)
mutationData$VAF_CN_corrected = mutationData$VAF

for(i in 1:nrow(mutationData)){
  if(mutationData$CN_status[i] == 1){
    mutationData$VAF_CN_corrected[i] = mutationData$VAF[i]/2
  }
  if(mutationData$CN_status[i] == 3){
    mutationData$VAF_CN_corrected[i] = mutationData$VAF[i]*1.5
  }
}

for(i in 1:nrow(mutationData)){
  if(!is.na(mutationData$VAF_CN_corrected[i])){
    if(mutationData$VAF_CN_corrected[i] > 100){
      mutationData$VAF_CN_corrected[i] = 100
    }
  }
}

# extract useful columns and append to mutation
mutationData = dplyr::select(mutationData, Study_Variant_ID, VAF_CN_corrected)

mutationData$Study_Variant_ID = as.numeric(mutationData$Study_Variant_ID)
joined_muts_clinical$Study_Variant_ID = as.numeric(joined_muts_clinical$Study_Variant_ID)

joined_muts_clinical = left_join(joined_muts_clinical, mutationData, by = "Study_Variant_ID")


# ======================= #
# Combine dataframes ####
# ======================= #
final_data_matrix_sub <- rbind.fill(combined2, joined_muts_clinical)
# make sure to make the VAF in the appropriate range
final_data_matrix_sub$VAF <- round(final_data_matrix_sub$VAF, 2)

# for some reason, SRSF2 used to be SFRS2 and was not changed in any of the datasets except BeatAML. Find the old annotations and change them here:
for(i in 1:nrow(final_data_matrix_sub)){
  if(final_data_matrix_sub$symbol[i] == "SFRS2"){
    final_data_matrix_sub$symbol[i] <- "SRSF2"
  }
}

colnames(final_data_matrix_sub)[1] <- "Sample"



# ==================================================================== #
# Add clinical annotations (mainly survival) to the mutation data ####
# ==================================================================== #

# TCGA
# clinical data
TCGA_survival <- read.table("~/Desktop/MetaAML_results/raw_data/data_clinical_patient.txt", header = T, stringsAsFactors = F, sep = "\t")

TCGA_survival2 <- dplyr::select(TCGA_survival, PATIENT_ID, OS_STATUS, OS_MONTHS, SEX, AGE, CYTOGENETIC_CODE_OTHER, RISK_MOLECULAR, BM_BLAST_PERCENTAGE, PB_BLAST_PERCENTAGE, WBC)
# need to convert tcga data from months to days
TCGA_survival2$OS_MONTHS <- (TCGA_survival2$OS_MONTHS*30)

# Beat AML (Tyner et al.)
BeatAML_survival2 <- dplyr::select(BeatAML_sample_data_types,  LabId, PatientId, vitalStatus, overallSurvival, consensus_sex, ageAtDiagnosis, Karyotype, ELN2017, `%.Blasts.in.BM`, `%.Blasts.in.PB`, WBC.Count, Hemoglobin, LDH, Platelet.Count)

# Stanford Majeti lab
# survival data
Stanford_survival <- read_excel("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/151102_Clinical Outcomes of All Patients without_PHI.xlsx")
# remove empty patient rows
Stanford_survival <- Stanford_survival[1:133,]
Stanford_survival2 <- dplyr::select(Stanford_survival, `Patient SU Number`, `Overall Survival (Event = Death)`, `Duration from diagnosis to D1`, Gender, `Age at Diagnosis`, Cytogenetics, 15)
names(Stanford_survival2)[1] = "Sample"

clin_annotations <- read_excel("~/Desktop/MetaAML_results/raw_data/Archive/160104_RM_AML_Sample_Info.xlsx")
# #
# # # dplyr::select columns of interest
labels <- dplyr::select(clin_annotations, Sample, `% Blasts`)
#
Stanford_survival2 = left_join(Stanford_survival2, labels, by = "Sample")


# Papaemmanuil et al.
Multistage_survival2 <- dplyr::select(Multistage_survival, PDID, Status, OS, gender, AOD, complex, M_Risk, BM_Blasts, PB_Blasts, wbc, HB, LDH, platelet)

# make all the headers the same
colnames(TCGA_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Age", "Cytogenetics", "Risk", "BM_blast_percent", "PB_blast_percent", "WBC")

colnames(BeatAML_survival2) <- c("Sample", "PatientId", "Censor", "Time_to_OS", "Sex", "Age", "Cytogenetics", "Risk", "BM_blast_percent", "PB_blast_percent", "WBC", "Hemoglobin", "LDH", "Platelet")
# BeatAML_survival2$Cohort <- "Tyner"

colnames(Stanford_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Age", "Cytogenetics", "Risk", "BM_blast_percent")
# Stanford_survival2$Cohort <- "Stanford"

colnames(Multistage_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Age", "Cytogenetics", "Risk", "BM_blast_percent", "PB_blast_percent", "WBC", "Hemoglobin", "LDH", "Platelet")

# Multistage_survival2$Cohort <- "Papaemmanuil"

mutation_matrix_survival <- rbind.fill(TCGA_survival2, BeatAML_survival2, Stanford_survival2, Multistage_survival2)

# homogenize the censor calls
for(i in 1:nrow(mutation_matrix_survival)){
  if(mutation_matrix_survival$Censor[i] == "LIVING" | mutation_matrix_survival$Censor[i] == "Alive"){
    mutation_matrix_survival$Censor[i] <- 0
  }
  if(mutation_matrix_survival$Censor[i] == "DECEASED" | mutation_matrix_survival$Censor[i] == "Dead"){
    mutation_matrix_survival$Censor[i] <- 1
  }
}

# remove cases where the time to OS is negative
# mutation_matrix_survival <- mutation_matrix_survival[mutation_matrix_survival$Time_to_OS >= 0, ]


# ======================================= #
# Combine mutation and clinical data ####
# ======================================= #
final_data_matrix <- left_join(final_data_matrix_sub, mutation_matrix_survival,  by = "Sample")

# homogenize vatiant type annotations
for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$symbol[i] == "FLT3" & final_data_matrix$variant_type[i] == "tandem_duplication"){
    final_data_matrix$variant_type[i] <- "ITD"
  }
  if(final_data_matrix$variant_type[i] == "deletion"){
    final_data_matrix$variant_type[i] <- "Deletion"
  }
  if(final_data_matrix$variant_type[i] == "insertion"){
    final_data_matrix$variant_type[i] <- "Insertion"
  }
}

final_data_matrix <- as.data.frame(final_data_matrix)

final_data_matrix$Risk[is.na(final_data_matrix$Risk)] <- "NA"

for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$Risk[i] == "Inter-1" | final_data_matrix$Risk[i] == "Inter-2" | final_data_matrix$Risk[i] == "FavorableOrIntermediate" | final_data_matrix$Risk[i] == "IntermediateOrAdverse" | final_data_matrix$Risk[i] == 1){
    final_data_matrix$Risk[i] <- "Intermediate"
  }
  if(final_data_matrix$Risk[i] == "Good" | final_data_matrix$Risk[i] == 0){
    final_data_matrix$Risk[i] <- "Favorable"
  }
  if(final_data_matrix$Risk[i] == "Poor" | final_data_matrix$Risk[i] == 2){
    final_data_matrix$Risk[i] <- "Adverse"
  }
  if(final_data_matrix$Risk[i] == "N.D." | final_data_matrix$Risk[i] == "NA"){
    final_data_matrix$Risk[i] <- "Unknown"
  }
}

final_data_matrix <- final_data_matrix[!duplicated(final_data_matrix), ]




# ========================== #
# Other smaller cohorts ####
# ========================== #

# For many of these cohorts, I manually scraped results from supplimental PDFs and converted them into excel files

# Blood_2014_Lindsley ####
lindsley1 <- read_excel("~/Desktop/MetaAML_results/Data/Blood_2014_Lindsley.xlsx")
lindsley2 <- read_excel("~/Desktop/MetaAML_results/Data/Blood_2014_Lindsley_additional_clinical.xlsx")

# format small details to merge data
lindsley1 <- dplyr::select(lindsley1, UPN, Gene, Result, VAF, AA)
colnames(lindsley1) <- c("Sample", "symbol", "variant_type", "VAF", "amino_acid_change")

lindsley2$SUBJID <- gsub("-", "_", lindsley2$SUBJID)
lindsley2 <- dplyr::select(lindsley2, 1:3,5,7:9)
colnames(lindsley2) <- c("Sample", "Age", "Sex", "Subset", "Censor", "Time_to_OS", "Cytogenetics")

lindsley <- left_join(lindsley1, lindsley2,  by = "Sample")

lindsley$Time_to_OS <- (lindsley$Time_to_OS*30)

lindsley$Cohort <- "Lindsley"

# annotate the type of mutation
for(i in 1:nrow(lindsley)){
  if(lindsley$symbol[i] == "FLT3" & lindsley$variant_type[i] == "nonsynonymous SNV"){
    lindsley$variant_type[i] = "SNV"
  }
  if(lindsley$symbol[i] == "FLT3" & lindsley$variant_type[i] == "nonframeshift deletion"){
    lindsley$variant_type[i] = "Deletion"
  }
  if(lindsley$symbol[i] == "FLT3" & lindsley$variant_type[i] == "nonframeshift insertion"){
    lindsley$variant_type[i] = "ITD"
  }
}

lindsley$VAF = lindsley$VAF*100

# n_distinct(lindsley$Sample)
# 183

# Oncotarget_2016_Wang ####
wang <- read_excel("~/Desktop/MetaAML_results/Data/Oncotarget_2016_Wang.xlsx")
wang <- wang[-1,]
wang <- dplyr::select(wang, `Sample ID`, Gene, `Amino acid`)
wang$Subset <- "de_novo"
colnames(wang)[1:4] <- c("Sample", "symbol", "amino_acid_change", "Subset")
wang$Cohort <- "Wang"
wang$variant_type = NA

for(i in 1:nrow(wang)){
  if(wang$symbol[i] == "FLT3"){
    wang$variant_type[i] = "SNV"
  }
  if(wang$symbol[i] == "FLT3-ITD"){
    wang$variant_type[i] = "ITD"
  }
}

# n_distinct(wang$Sample)
# 91

# BMC_2016_Au ####
Au <- read_excel("~/Desktop/MetaAML_results/Data/BMC_2016_Au.xlsx")
Au <- Au %>% separate(`Sex/Age`, sep = "/", into = c("Sex", "Age"))
Au <- Au[,c(1:4,6,8,9)]
colnames(Au) <- c("Sample", "Sex", "Age", "Subset", "Cytogenetics", "symbol", "amino_acid_change")
Au$variant_type = NA

Au$Cohort <- "Au"
# n_distinct(Au$Sample)
# 50


# NEJM_2016_Welch ####
welch_mut <- read_excel("~/Desktop/MetaAML_results/Data/NEJM_2016_Welch_appendix_variants.xlsx")
welch_mut <- welch_mut[,c(1:3, 8, 16)]
colnames(welch_mut)[1:5] <- c("Sample", "symbol", "amino_acid_change", "variant_type", "VAF")
# welch_mut$VAF <- (welch_mut$VAF/100)

welch_clinical <- read_excel("~/Desktop/MetaAML_results/Data/NEJM_2016_Welch_clinical.xlsx")
welch_clinical <- na.omit(dplyr::select(welch_clinical, 1,3,5:9,18:21))
colnames(welch_clinical) <- c("Sample", "Subset", "Age", "Sex", "PB_wbc_percent", "PB_blast_percent", "BM_blast_percent", "Time_to_OS", "Censor", "Risk", "Cytogenetics")

# combind mutations and clinical together
welch <- left_join(welch_mut, welch_clinical, by = "Sample")

for(i in 1:nrow(welch)){
  if(!is.na(welch$Risk[i])){
    if(welch$Risk[i] == 0){
      welch$Risk[i] <- "Favorable"
    }
    if(welch$Risk[i] == 1){
      welch$Risk[i] <- "Intermediate"
    }
    if(welch$Risk[i] == 2){
      welch$Risk[i] <- "Adverse"
    }
  }
}

welch$Cohort <- "Welch"

# n_distinct(welch$Sample)
# 56
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

final_data_matrix <- rbind.fill(final_data_matrix, lindsley, wang, Au, welch)


final_data_matrix$variant_type[is.na(final_data_matrix$variant_type)] <- "Unknown"

# make all the variant type annotations match
for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$variant_type[i] == "frameshift insertion" | final_data_matrix$variant_type[i] == "nonframeshift insertion" | final_data_matrix$variant_type[i] == "INS"){
    final_data_matrix$variant_type[i] <- "Insertion"
  }
  if(final_data_matrix$variant_type[i] == "stopgain SNV" | final_data_matrix$variant_type[i] == "nonsynonymous SNV" | final_data_matrix$variant_type[i] == "SNP"){
    final_data_matrix$variant_type[i] <- "SNV"
  }
  if(final_data_matrix$variant_type[i] == "DEL" | final_data_matrix$variant_type[i] == "frameshift deletion" | final_data_matrix$variant_type[i] == "nonframeshift deletion"){
    final_data_matrix$variant_type[i] <- "Deletion"
  }
  if(final_data_matrix$variant_type[i] == "exonic;splicing" | final_data_matrix$variant_type[i] == "splicing"){
    final_data_matrix$variant_type[i] <- "Splicing"
  }
}

final_data_matrix$VAF <- as.numeric(final_data_matrix$VAF)

# correct any FLT3 ITD differences in annotation
for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$symbol[i] == "FLT3-ITD"){
    final_data_matrix$symbol[i] <- "FLT3"
  }
}



#  Blood_2015_Garg ####
# downloaded the supplimental tables and manually removed false mutations and germline mutations to make the final file
garg_clin_1 <- read_excel("~/Desktop/MetaAML_results/Data/Garg_2015_blood_clinical_discovery_cohort.xlsx")
garg_mut_1 <- read_excel("~/Desktop/MetaAML_results/Data/Garg_2015_blood_mutations_discovery_cohort.xlsx")
garg_clin_2 <- read_excel("~/Desktop/MetaAML_results/Data/Garg_2015_blood_clinical_targeted_cohort.xlsx")
garg_mut_2 <- read_excel("~/Desktop/MetaAML_results/Data/Garg_2015_blood_mutations_targeted_cohort.xlsx")

colnames(garg_clin_2) <- colnames(garg_clin_1)
garg_clin <- rbind(garg_clin_1, garg_clin_2)

garg_mut <- rbind(garg_mut_1, garg_mut_2)

garg_data <- left_join(garg_mut, garg_clin, by = "Sample ID")

# converty the OS data to days
garg_data[,20] <- as.numeric(unlist(garg_data[,20]))
garg_data[,20] <- garg_data[,20]*30

# dplyr::select only diagnostic mutations
garg_data <- garg_data[grep("DX", garg_data$`Status in Disease`), ]
# split the age/sex column
garg_data <- garg_data %>% separate(`Age/ Sex`, sep = "/", into = c("Age", "Sex"))

# subset to useful columns
garg_data <- dplyr::select(garg_data, `Sample ID`, `Annotated Gene`, `Mutation Type`, `Amino Acid Change`, 21, `Dead/ Alive`, Sex, Age, `Cytogenetic (DX)`)

colnames(garg_data) <- c("Sample", "symbol", "variant_type", "amino_acid_change", "Time_to_OS", "Censor", "Sex", "Age", "Cytogenetics")

garg_data$Subset <- "de_novo"

# split the age/sex column
flt3_add <- garg_clin %>% separate(`Age/ Sex`, sep = "/", into = c("Age", "Sex"))

# subset to useful columns
flt3_add <- dplyr::select(flt3_add, `Sample ID`, `Dead/ Alive`, Sex, Age, `Cytogenetic (DX)`, 14)
flt3_add$symbol <- "FLT3"
# manually add the ITD calls (since all patients have ITD in this cohort)
flt3_add$variant_type <- "ITD"
flt3_add$Subset <- "de_novo"

colnames(flt3_add)[1] <- "Sample"
colnames(flt3_add)[2] <-"Censor"
colnames(flt3_add)[5] <-"Cytogenetics"
colnames(flt3_add)[6] <-"Time_to_OS"

flt3_add$Time_to_OS <- as.numeric(flt3_add$Time_to_OS)

flt3_add$Time_to_OS <- unlist(flt3_add$Time_to_OS)*30

# add the flt3 mutations
garg_data <- rbind.fill(garg_data, flt3_add)

# homogenize the annotations with the other datasets
for(i in 1:nrow(garg_data)){
  if(garg_data$Censor[i] == "Alive"){
    garg_data$Censor[i] <- 0
  }
  if(garg_data$Censor[i] == "Dead"){
    garg_data$Censor[i] <- 1
  }
  if(garg_data$Sex[i] == "M"){
    garg_data$Sex[i] <- "Male"
  }
  if(garg_data$Sex[i] == "F"){
    garg_data$Sex[i] <- "Female"
  }
  if(garg_data$variant_type[i] == "Splice site"){
    garg_data$variant_type[i] <- "Splicing"
  }
  if(garg_data$variant_type[i] == "Frameshift deletion" | garg_data$variant_type[i] == "Non Frameshift deletion"){
    garg_data$variant_type[i] <- "Deletion"
  }
  if(garg_data$variant_type[i] == "Frameshift\r\ninsertion" | garg_data$variant_type[i] == "Frameshift\r\ndeletion" | garg_data$variant_type[i] == "Frameshift insertion"){
    garg_data$variant_type[i] <- "INDEL"
  }
  if(garg_data$variant_type[i] == "Stoploss" | garg_data$variant_type[i] == "stopgain" | garg_data$variant_type[i] == "Stopgain" | garg_data$variant_type[i] == "Missense" ){
    garg_data$variant_type[i] <- "SNV"
  }
}

garg_data$Cohort <- "Garg"

final_data_matrix <- rbind.fill(final_data_matrix, garg_data)


# Huet et al Blood 2018 ####
# download supplimental PDFs and manually extract data and put into an excell file
download.file("https://ash.silverchair-cdn.com/ash/content_public/journal/blood/132/8/10.1182_blood-2018-03-840348/4/blood840348-sup-tables2.xlsx?Expires=1579153943&Signature=VCqH4IQOtZ~x9YPe1dJSz1AHegirMKv6N4euRcbvd3gsayA~RULg17B6mjZY0lrUzbbbj5D7MdLJNzI0sUk4QFDz6IKvvWsV~l1CQL4J1lOHejvckLQzXikTy-4SD62jwmDMgSeBMhXqVo6pmdHXI0iO7b8ARJV8EDZ9Bt8h~h6nIXOaf5JZqe8hehacxWIMjCWRd5OulQjpqQu10MLfK6GrpISgRGDYe1WODtkPIQ~o6qQj~eqTFGvzwbSmhWuRGSAFKEFzPzgYLiKmkQ8errjq9vI2uHgpurlwtB98gQOnizt9goXVoDPsGd0P~fk8CPPEGJdUtSkOB41b8p3JTg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA", destfile = "~/Desktop/MetaAML_results/raw_data/blood840348-sup-tables2.xlsx")

huet=read_excel("~/Desktop/MetaAML_results/Data/Huet_2018_Blood.xlsx")
colnames(huet)=huet[1,]
huet=huet[-1,]

# transform the data for easier appending to MetaAML dataframe
muts=huet[c(1,41:99)]
colnames(muts)[56]="U2AF1"

for(i in 1:nrow(muts)){
  muts$`Patient ID`[i] = paste("pt", i, "huet", sep = "_")
}
muts = data.table(muts)

muts_m=melt(data = muts, id.vars = "Patient ID")

muts_m$variable=gsub(" .*", "", muts_m$variable)
muts_m=subset(muts_m, muts_m$value == 1)

muts_m=muts_m[,c(1:2)]
colnames(muts_m)[2]="symbol"

muts_m$symbol = sub("\\?_","\\?",muts_m$symbol)

# subset clinical data
clin=huet[c(1:13,100:103)]

# add unique patient ID to ensure sample uniqueness
for(i in 1:nrow(clin)){
  clin$`Patient ID`[i] = paste("pt", i, "huet", sep = "_")
}

# merge data together
huet=left_join(clin, muts_m, by = "Patient ID")
huet$Cohort="Huet"

# clean up column and cell annotations to match the MetaAML format
# column names
setnames(huet, c("Patient ID", "de novo, secondary or therapy related AML", "ELN 2017 category", "Absolute leucocyte count (109/L)", "LDH (UI/L)", "Hemoglobin (gr/L)", "Platelet count (109/L)", "Blast cells in bone marrow (%)", "Blast cells in peripheral blood (%)", "Time to final outcome (months)", "Final outcome"), c("Sample", "Subset", "Risk", "WBC", "LDH", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent", "Time_to_OS", "Censor"))

# cells
huet$Sex <- ifelse(huet$Sex == "F", "Female","Male")

for(i in 1:nrow(huet)){
  if(huet$Risk[i] == "Fav"){
    huet$Risk[i] = "Favorable"
  }
  if(huet$Risk[i] == "Int"){
    huet$Risk[i] = "Intermediate"
  }
  if(huet$Risk[i] == "Adv"){
    huet$Risk[i] = "Adverse"
  }
  if(huet$Censor[i] %in% c("Alive in first CR","Alive after relapse")){
    huet$Censor[i] = 0
  }
  if(huet$Censor[i] %in% c("Non remission death","Relapse death", "Non relapse death")){
    huet$Censor[i] = 1
  }
}

huet = huet %>% separate(symbol, c("symbol", "amino_acid_change"), sep = "-")
huet = huet %>% separate(symbol, c("symbol", "type_2"), sep = "_")

for(i in 1:nrow(huet)){
  if(!is.na(huet$type_2[i])){
    huet$amino_acid_change[i]=huet$type_2[i]
  }
}
huet$type_2=NULL

huet$Time_to_OS=as.numeric(huet$Time_to_OS)

huet$Time_to_OS <- huet$Time_to_OS * 30

huet$variant_type = NA

huet$symbol = as.character(huet$symbol)

for(i in 1:nrow(huet)){
  if(!is.na(huet$symbol[i])){
    if(huet$symbol[i] == "FLT3"){
      if(!is.na(huet$amino_acid_change[i])){
        huet$variant_type[i] = huet$amino_acid_change[i]
      }
    }
  }
}

final_data_matrix <- rbind.fill(final_data_matrix, huet)


# Hirsch et al 2016 Nat Com ####
# supplimental data pdf files were downloaded, converted to excel files, and manually formatted
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12475/MediaObjects/41467_2016_BFncomms12475_MOESM1148_ESM.xlsx", destfile = "~/Desktop/MetaAML_results/raw_data/41467_2016_BFncomms12475_MOESM1148_ESM.xlsx")
#
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fncomms12475/MediaObjects/41467_2016_BFncomms12475_MOESM1148_ESM.xlsx", destfile = "~/Desktop/MetaAML_results/raw_data/41467_2016_BFncomms12475_MOESM1149_ESM.xlsx")

hirsch_clinical=read_excel("~/Desktop/MetaAML_results/Data/Hirsch_2016_Nat_Comm_clinical_data.xlsx", sheet = 1)
hirsch_muts_diagnosis=read_excel("~/Desktop/MetaAML_results/Data/Hirsch_2016_Nat_Comm_mutations.xlsx", sheet = 1)
hirsch_muts_relapse=read_excel("~/Desktop/MetaAML_results/Data/Hirsch_2016_Nat_Comm_mutations.xlsx", sheet = 2)

# format data for apending to maset dataframe
setnames(hirsch_clinical, c("UPN", "Gender", "WBC (x 10 9/L)", "% Blast cells in BM", "time from diagnosis (days)", "alive at last follow-up"), c("Sample", "Sex", "WBC", "BM_blast_percent", "Time_to_OS", "Censor"))

hirsch_clinical$Censor <- ifelse(hirsch_clinical$Censor == "yes", 0,1)

hirsch_clinical$Subset = NA

# assign type of AML annotation to match the MetaAML format
for(i in 1:nrow(hirsch_clinical)){
  if(hirsch_clinical$WHO[i] %in% c("AML NOS","AML with NPM1 mutation","AML with NPM1\r\nmutation","AML with t(8;21)","AML with inv (16)")){
    hirsch_clinical$Subset[i] = "de_novo"
  }
  if(hirsch_clinical$WHO[i] == "sAML"){
    hirsch_clinical$Subset[i] = "secondary"
  }
  if(hirsch_clinical$WHO[i] == "tAML"){
    hirsch_clinical$Subset[i] = "therapy"
  }
}

# add mutation data to clinical info
hirsch_muts_relapse$diagnosis_or_relapse="relapse"
hirsch_muts_diagnosis$diagnosis_or_relapse="diagnosis"
colnames(hirsch_muts_relapse)=colnames(hirsch_muts_diagnosis)
hirsch_muts_all=rbind(hirsch_muts_diagnosis, hirsch_muts_relapse)

colnames(hirsch_muts_all)=c("Sample", "symbol", "Chromosome", "location", "ref base", "alt base", "VAF", "amino_acid_change", "NM", "diagnosis_or_relapse")

hirsch_all=left_join(hirsch_muts_all, hirsch_clinical, by = "Sample")

hirsch_all$variant_type = "Unknown"

for(i in 1:nrow(hirsch_all)){
  if(hirsch_all$symbol[i] == "FLT3"){
    if(nchar(hirsch_all$`ref base`)[i] > 1 | nchar(hirsch_all$`alt base`)[i] > 1){
      hirsch_all$variant_type[i] = "ITD"
    }
    if(nchar(hirsch_all$`ref base`)[i] == 1 & nchar(hirsch_all$`alt base`)[i] == 1){
      hirsch_all$variant_type[i] = "SNV"
    }
  }
}

hirsch_all$VAF = hirsch_all$VAF*100

hirsch_all$Cohort="Hirsch"

final_data_matrix <- rbind.fill(final_data_matrix, hirsch_all)


# Azizi ####
# Armon Azizi has curated a paired diagnostic and relapse cohort
clin_dat <- read_excel("~/Desktop/MetaAML_results/Data/azizi_diagnosis_relapse_aml_clinical.xlsx")
mut_dat <- read_excel("~/Desktop/MetaAML_results/Data/azizi_diagnosis_relapse_all_somatic_mutations.xlsx")

azizi_data <- left_join(mut_dat, clin_dat, by = "Case")

azizi_data <- dplyr::select(azizi_data, Case, Gene, DX_VAF, REL_VAF, AA_Change, Overall_Survival, Gender, Age, DX_Karyotype, ELN, AML_Type, DX_Blasts, REL_Blasts, Cohort)
colnames(azizi_data) <- c("Sample", "symbol", "Diagnostic_VAF","Relapse_VAF", "amino_acid_change", "Time_to_OS", "Sex", "Age", "DX_Cytogenetics", "Risk", "Subset", "DX_Blasts", "REL_Blasts", "Cohort")

# remove BeatAML data that is redundant from this cohort
azizi_data <- subset(azizi_data, Cohort!="beatAML")

# correct any FLT3 ITD differences in annotation
for(i in 1:nrow(azizi_data)){
  if(azizi_data$symbol[i] == "FLT3" & azizi_data$amino_acid_change[i] == "ITD"){
    azizi_data$symbol[i] <- "FLT3-ITD"
  }
  if(azizi_data$Diagnostic_VAF[i] == "P"){
    azizi_data$Diagnostic_VAF[i] <- NA
  }
  if(azizi_data$Relapse_VAF[i] == "P"){
    azizi_data$Relapse_VAF[i] <- NA
  }
}

# make duplicate dataframes for the diagnostic and relapse cases/vafs then append them together
azizi_diagnosis <- azizi_data[,c(1:3,5:12,14)]
colnames(azizi_diagnosis)[c(3,8,11)] <- c("VAF", "Cytogenetics", "BM_blast_percent")

# remove rows where the diagnosis vaf is 0
for(i in 1:nrow(azizi_diagnosis)){
  if(!is.na(azizi_diagnosis$VAF[i])){
    if(azizi_diagnosis$VAF[i] == 0){
      azizi_diagnosis <- azizi_diagnosis[-i,]
    }
  }
}

azizi_relapse <- azizi_data[,c(1:2,4:8,10,12,14)]
colnames(azizi_relapse)[c(3,9)] <- c("VAF", "BM_blast_percent")
azizi_relapse$Subset <- "relapse"

# remove rows where the relapse vaf is 0
for(i in 1:nrow(azizi_relapse)){
  if(!is.na(azizi_relapse$VAF[i])){
    if(azizi_relapse$VAF[i] == 0){
      azizi_relapse <- azizi_relapse[-i,]
    }
  }
}

azizi_data <- rbind.fill(azizi_diagnosis, azizi_relapse)

# make sure the NPM1 aa change annotations are INDELS
for(i in 1:nrow(azizi_data)){
  if(is.na(azizi_data$amino_acid_change[i]) &  azizi_data$symbol[i] == "NPM1"){
    azizi_data$amino_acid_change[i] <- "INDEL"
  }
}

# annotate the types of FLT3 mutations
azizi_data$variant_type = "Unknown"

for(i in 1:nrow(azizi_data)){
  if(azizi_data$symbol[i] == "FLT3" & azizi_data$amino_acid_change[i] == "ITD"){
    azizi_data$variant_type[i] = "ITD"
  }
  if(azizi_data$symbol[i] == "FLT3" & nchar(azizi_data$amino_acid_change)[i] > 10 ){
    azizi_data$variant_type[i] = "ITD"
  }
  if(azizi_data$symbol[i] == "FLT3" & nchar(azizi_data$amino_acid_change)[i] > 6 & nchar(azizi_data$amino_acid_change)[i] < 9){
    azizi_data$variant_type[i] = "SNV"
  }
}

azizi_data$DX_VAF = azizi_data$DX_VAF*100

azizi_data$Cohort = azizi_data$Cohort[azizi_data$Cohort == "Grief"] <- "Greif"

final_data_matrix <- rbind.fill(final_data_matrix, azizi_data)

# homogenize data coding and annotations ####

# homogenize the annnotations for type of AML
for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$Subset[i])){
    if(final_data_matrix$Subset[i] == "tAML" | final_data_matrix$Subset[i] == "TAML"| final_data_matrix$Subset[i] == "therapy related"){
      final_data_matrix$Subset[i] <- "therapy"
    }
    if(final_data_matrix$Subset[i] == "transformed" | final_data_matrix$Subset[i] == "sAML"){
      final_data_matrix$Subset[i] <- "secondary"
    }
    if(final_data_matrix$Subset[i] == "AML" | final_data_matrix$Subset[i] == "de novo"){
      final_data_matrix$Subset[i] <- "de_novo"
    }
    if(final_data_matrix$Subset[i] == "Relapsed_AML"){
      final_data_matrix$Subset[i] <- "relapse"
    }
  }
}

final_data_matrix$variant_type[is.na(final_data_matrix$variant_type)] <- "Unknown"
final_data_matrix$VAF <- as.numeric(final_data_matrix$VAF)
final_data_matrix$Time_to_OS <- as.numeric(final_data_matrix$Time_to_OS)

# sex coding
for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$Sex[i])){
    if(final_data_matrix$Sex[i] == "M" | final_data_matrix$Sex[i] == 1){
      final_data_matrix$Sex[i] <- "Male"
    }
    if(final_data_matrix$Sex[i] == "F" | final_data_matrix$Sex[i] == 2){
      final_data_matrix$Sex[i] <- "Female"
    }
  }
}

# correct variant type annotations manually
for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$amino_acid_change[i])){
    if(final_data_matrix$symbol[i] == "FLT3" & final_data_matrix$amino_acid_change[i] == "ITD"){
      final_data_matrix$variant_type[i] <- "ITD"
    }
    if(final_data_matrix$symbol[i] == "FLT3" & final_data_matrix$variant_type[i] == "Insertion"){
      final_data_matrix$variant_type[i] <- "ITD"
    }
    if(final_data_matrix$amino_acid_change[i] == "FS"){
      final_data_matrix$variant_type[i] <- "INDEL"
    }
    if(final_data_matrix$variant_type[i] == "TKD"){
      final_data_matrix$variant_type[i] <- "SNV"
    }
  }
}

# remove CML patients
final_data_matrix <- final_data_matrix[!grepl("CML", final_data_matrix$Subset),]

# remove SBTB33 mutations because they are over-represented in the Majeti cohort
final_data_matrix = subset(final_data_matrix, final_data_matrix$symbol != "ZBTB33")

# frequency of gene mutations
pt_gene <- as.data.frame(unique(final_data_matrix[,c(1:2)]))
pt_gene$mut_freq_gene <- NA

for(i in 1:nrow(pt_gene)){
  pt_gene_sub <- subset(pt_gene, pt_gene$symbol == pt_gene[i,2])
  pt_gene$mut_freq_gene[i] <- as.numeric(nrow(pt_gene_sub))
}

# subset final matrix to only include mutations called in at least two patients
pt_gene <- subset(pt_gene, pt_gene$mut_freq_gene > 2)

final_data_matrix <- setDT(final_data_matrix)[symbol %chin% pt_gene$symbol]
final_data_matrix <- final_data_matrix %>% left_join(pt_gene, by=c("Sample","symbol"))

# finally, add a column with the number of unique mutations per patient
# loop through the meta file and count the number of mutations per patient and append a column with the number per patient

pt_mut <- as.data.frame(unique(final_data_matrix[,c(1:2)]))
pt_mut$mut_freq_pt <- NA

for(i in 1:nrow(pt_mut)){
  pt_gene_sub <- subset(pt_mut, pt_mut$Sample == pt_mut[i,1])
  pt_mut$mut_freq_pt[i] <- as.numeric(nrow(pt_gene_sub))
}

final_data_matrix <- final_data_matrix %>% left_join(pt_mut, by=c("Sample","symbol"))

# final_data_matrix <-  final_data_matrix[!duplicated(final_data_matrix[1:5]),]
# add a column for the mutation bin per patient
final_data_matrix$mut_freq_bin = NA

for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$mut_freq_pt[i] == 1){
    final_data_matrix$mut_freq_bin[i] <- 1
  }
  if(final_data_matrix$mut_freq_pt[i] == 2){
    final_data_matrix$mut_freq_bin[i] <- 2
  }
  if(final_data_matrix$mut_freq_pt[i] == 3){
    final_data_matrix$mut_freq_bin[i] <- 2
  }
  if(final_data_matrix$mut_freq_pt[i] == 3){
    final_data_matrix$mut_freq_bin[i] <- 3
  }
  if(final_data_matrix$mut_freq_pt[i] == 4){
    final_data_matrix$mut_freq_bin[i] <- 4
  }
  if(final_data_matrix$mut_freq_pt[i] == 5){
    final_data_matrix$mut_freq_bin[i] <- 5
  }
  if(final_data_matrix$mut_freq_pt[i] == 6){
    final_data_matrix$mut_freq_bin[i] <- 6
  }
  if(final_data_matrix$mut_freq_pt[i] == 7){
    final_data_matrix$mut_freq_bin[i] <- 7
  }
  if(final_data_matrix$mut_freq_pt[i] > 7){
    final_data_matrix$mut_freq_bin[i] <- "8+"
  }
}

# # correct VAFs based on x-linkage and gender
# download a file of X-linked genes
# https://www.ncbi.nlm.nih.gov/gene/?term=X%5BCHR%5D+AND+human%5BORGN%5D
x_genes <- read.table("~/Downloads/gene_result.txt", sep = "\t", header = T, fill = T,  stringsAsFactors = FALSE, quote = "")
x_genes$Symbol <- as.character(x_genes$Symbol)

final_data_matrix$VAF_male_x <- NA

`%!in%` = Negate(`%in%`)

for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$Sex[i])){
    if(!is.na(final_data_matrix$VAF[i])){
      if(final_data_matrix$Sex[i] == "Male"){
        if(final_data_matrix$symbol[i] %chin% x_genes$Symbol){
          final_data_matrix$VAF_male_x[i] <- (final_data_matrix$VAF[i])/2
        } else {
          final_data_matrix$VAF_male_x[i] <- final_data_matrix$VAF[i]
        }
      }
      if(final_data_matrix$Sex[i] == "Female"){
        final_data_matrix$VAF_male_x[i] <- final_data_matrix$VAF[i]
      }
    }
  }
}

colnames(final_data_matrix)[2] <- "Gene"
colnames(final_data_matrix)[7] <- "PatientId"

final_data_matrix$PatientId.y <- NULL
# final_data_matrix[,10:13] <- NULL

final_data_matrix <- final_data_matrix %>% distinct(Sample, Gene, amino_acid_change, .keep_all = TRUE)

for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "Insertion"){
    final_data_matrix$variant_type[i] <- "ITD"
  }
}

for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$Censor[i])){
    if(final_data_matrix$Censor[i] == "Unknown"){
      final_data_matrix$Censor[i] <- NA
    }
  }
}


# finalize the risk annotations
for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$Risk[i])){
    if(final_data_matrix$Risk[i] == "Intermediate-I"){
      final_data_matrix$Risk[i] <- "Intermediate"
    }
    if(final_data_matrix$Risk[i] == "Unfavorable"){
      final_data_matrix$Risk[i] <- "Adverse"
    }
  }
}

# complex karyotype
for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$Cohort[i] == "Papaemmanuil" & final_data_matrix$Cytogenetics[i] == 1){
    final_data_matrix$Cytogenetics[i] = "Complex Cytogenetics"
  }
  if(final_data_matrix$Cohort[i] == "Papaemmanuil" & final_data_matrix$Cytogenetics[i] == 0){
    final_data_matrix$Cytogenetics[i] = "Normal Karyotype"
  }
}

# Make the stanford data Majeti lab annotated
for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$Cohort[i] == "Stanford"){
    final_data_matrix$Cohort[i] = "Majeti"
  }
}

# add a column for annotating the mutation category
final_data_matrix$mutation_category <- NA

DNA_methylation <- list("DNMT3A","IDH2","TET2","IDH1")
Chromatin_cohesin <- list("ASXL1", "RAD21", "STAG2", "EZH2", "BCOR")
RTK_RAS_Signaling <- list("PTPN11", "CBL", "NF1", "KRAS", "KIT", "NRAS", "FLT3-ITD", "FLT3-TKD")
Splicing <- list("SF3B1", "SRSF2", "U2AF1")
Transcription <- list("CEBPA", "GATA2", "RUNX1", "MYC", "ETV6", "ZBTB33")
Tumor_suppressors <- list("TP53", "PHF6", "WT1")

for(i in 1:nrow(final_data_matrix)){
  if(final_data_matrix$Gene[i] %in% DNA_methylation){
    final_data_matrix$mutation_category[i] <- "DNA Methylation"
  }
  if(final_data_matrix$Gene[i] %in% Chromatin_cohesin){
    final_data_matrix$mutation_category[i] <- "Chromatin/Cohesin"
  }
  if(final_data_matrix$Gene[i] %in% RTK_RAS_Signaling){
    final_data_matrix$mutation_category[i] <- "RTK/RAS Signaling"
  }
  if(final_data_matrix$Gene[i] %in% Splicing){
    final_data_matrix$mutation_category[i] <- "Splicing"
  }
  if(final_data_matrix$Gene[i] %in% Transcription){
    final_data_matrix$mutation_category[i] <- "Transcription"
  }
  if(final_data_matrix$Gene[i] == "NPM1"){
    final_data_matrix$mutation_category[i] <- "NPM1"
  }
  if(final_data_matrix$Gene[i] %in% Tumor_suppressors){
    final_data_matrix$mutation_category[i] <- "Tumor suppressors"
  }
}

# dplyr::select final columns and save the data frame

final_data_matrix =  final_data_matrix %>%
  dplyr::select(Sample,Gene,VAF,variant_type,amino_acid_change,Subset,PatientId,specimenType,Cohort,Censor,Time_to_OS,Sex,Age,Cytogenetics,Risk,BM_blast_percent,PB_blast_percent,WBC,Hemoglobin,LDH,Platelet,PB_wbc_percent,mut_freq_gene,mut_freq_pt,mut_freq_bin,VAF_male_x,VAF_CN_corrected,mutation_category)

save(final_data_matrix,  file = "~/Desktop/MetaAML_results/final_data_matrix.RData")
write.csv(final_data_matrix,  file = "~/Desktop/MetaAML_results/final_data_matrix.csv")

# clean up the user's environment from intermediate files
rm(list=setdiff(ls(), "final_data_matrix"))




# ========================================================================================================================================== #
# Figure 1A ####
# UpSet plot of cohort ####
# MetaAML_data_upset_plot
# 
# This script plots the overlap of different data types for the meta analysis cohort
#
# load required packages
if (!require('UpSetR')) install.packages('UpSetR'); library('UpSetR')

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


p = upset(upset_matrix,
          nsets = 7,
          order.by = "freq",
          # main.bar.color = "#8073ac",
          # sets.bar.color = "",
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

# This script summarizes the major genetic features of the Meta-AML

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
# if (!require('ComplexHeatmap')) install.packages('ComplexHeatmap'); library('ComplexHeatmap')
library(devtools)
install_github("jokergoo/ComplexHeatmap")
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

# plot only genes mutated at least 100 times in the cohort
final_data_matrix_3 <- subset(final_data_matrix_3, final_data_matrix_3$mut_freq_gene >= 100)

# c = as.numeric(n_distinct(final_data_matrix_3$Sample))
c = 2865
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

png(filename = "~/Desktop/MetaAML_results/Figure_1/cohort_oncoprint_100.png", res = 300, width = 15, height = 6.5, units = "in")

fig1B
print(fig1B)
dev.off()
