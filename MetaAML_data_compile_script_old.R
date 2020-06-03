# metaAML_data_sompile_script server.R
#
# Brooks Benard
# bbenard@stanford.edu
#
#

### This script compiles the mutation calls across four different cohorts of AML patients (Stanford, TCGA, BeatAML, and MultiStage) and analyses the co-occurence, mutual exclusive, and survival patterns for different mutation pairs. It allows the user to perform subset analyses on de novo, secondary, and repalse patient cohorts, in addition to subsetting analysis to desired gene sets.

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
if (!require('survminer')) install.packages('survminer'); library('survminer')
#
# make directory dynamic per user
dir.create("~/Desktop/MetaAML_results/raw_data")
setwd("~/Desktop/MetaAML_results/raw_data/")
#
# download data from BeatAML website
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx")
#
#
# download data from TCGA website
# mutation data
download.file("http://download.cbioportal.org/laml_tcga_pub.tar.gz", destfile = "~/Desktop/MetaAML_results/raw_data/laml_tcga_pub.tar.gz")
untar("~/Desktop/MetaAML_results/raw_data/laml_tcga_pub.tar.gz")
TCGA_survival <- read.table("~/Desktop/MetaAML_results/raw_data/data_clinical_patient.txt", header = T, stringsAsFactors = F, sep = "\t")
#
download.file("https://tcga-data.nci.nih.gov/docs/publications/laml_2012/SupplementalTable06.tsv", destfile = "~/Desktop/MetaAML_results/raw_data/supplementalTable06.tsv")
#
#
# clinical annotations for Multistage
download.file("https://github.com/gerstung-lab/AML-multistage/blob/master/data/AMLSG_Clinical_Anon.RData?raw=true", destfile = "~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")
load("~/Desktop/MetaAML_results/raw_data/AMLSG_Clinical_Anon.RData")
write.table(clinicalData, "~/Desktop/MetaAML_results/raw_data/AML_knowledge_bank_data_clinical.txt")
#
#
# download the Multistage mutations
x <- getURL("https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/AMLSG_Genetic.txt")
#
#
# Stanford survival data
Stanford_survival <- read_excel("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/151102_Clinical Outcomes of All Patients without_PHI.xlsx")
#
#
#
# Let et al data
download.file("https://www.cell.com/cms/10.1016/j.cell.2012.06.023/attachment/8e4ea03d-9fb3-4701-8700-fce5e4ea96a8/mmc3.xls", destfile = "~/Desktop/MetaAML_results/raw_data/mmc3.xls")
ley_mut <- read_excel("~/Desktop/MetaAML_results/raw_data/mmc3.xls")
colnames(ley_mut) = ley_mut[1, ] # the first row will be the header
ley_mut = ley_mut[-1, ]   
# clinical data
download.file("https://www.cell.com/cms/10.1016/j.cell.2012.06.023/attachment/aa04d227-b485-42d3-9cdb-35615e02c45d/mmc1.xls", destfile = "~/Desktop/MetaAML_results/raw_data/mmc1.xls")
ley_clin <- read_excel("~/Desktop/MetaAML_results/raw_data/mmc1.xls")
colnames(ley_clin) = ley_clin[1, ] # the first row will be the header
ley_clin = ley_clin[-1, ]   

#
#
#
#


## Using the different datasets, create a single matrix containing all patients, variables, and values possible. This data frame will be used as the input for the various functions used to plot the results

#### load, subset, and arrange all the datasets for final compilation ####
#### TCGA: mutation matrix ####
# read in TCGA variants for analysis data
TCGAL_variants_raw <-read.table("~/Desktop/MetaAML_results/raw_data/supplementalTable06.tsv", header = T, sep = "\t", stringsAsFactors = F)

# subset to useful columns
variants <- TCGAL_variants_raw %>%
  select("TCGA_id", "gene_name", "type", "trv_type", "TumorVAF", "amino_acid_change")

# remove silent mutations and rows containing "-"
variants <- variants[!grepl("silent", variants$trv_type),]
variants <- variants[!grepl("intronic", variants$trv_type),]


# read in TCGA variants for analysis data
TCGAL_variants <-read.table("~/Desktop/MetaAML_results/raw_data/data_mutations_extended.txt", header = T, sep = "\t", stringsAsFactors = F)

# subset to useful columns
variants_2 <- TCGAL_variants %>%
  select("Hugo_Symbol", "Variant_Type")

colnames(variants_2)[2] <- "variant_type"

# count the number of mutations per gene and filter to recurrent genes
mut_table_2 <- aggregate(data.frame(count = variants_2), list(value = variants_2$Hugo_Symbol), length)
mut_table_2 <- select(mut_table_2, "value", "count.Hugo_Symbol")
mut_table_2 <- subset(mut_table_2, mut_table_2$count.Hugo_Symbol > 2)

# subset to recurrent mutations
variants <- setDT(variants)[gene_name %chin% mut_table_2$value] 

mut_table <- aggregate(data.frame(count = variants), list(value = variants$gene_name), length)
mut_table <- select(mut_table, "value", "count.gene_name")
colnames(mut_table)[1] <- "gene_name"

variants <- dplyr::left_join(variants, mut_table, by = "gene_name")

# count the number of mutations per patient and add to recurrent genes
mut_table_pts <- aggregate(data.frame(count = variants), list(value = variants$TCGA_id), length)
mut_table_pts <- select(mut_table_pts, "value", "count.TCGA_id")
colnames(mut_table_pts)[1] <- "TCGA_id"

# combine for final data frame
mut_table_final_TCGA <- dplyr::left_join(variants, mut_table_pts, by = "TCGA_id")
mut_table_final_TCGA$TumorVAF <- (mut_table_final_TCGA$TumorVAF/100)

# make column identifying subset type (i.e. de novo)
mut_table_final_TCGA$Subset <- "de_novo"
#
#
#
#
#
#
#

#### BeatAML: prep samples ####
# read in sample types
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# read in BeatAML variants for analysis data
BeatAML_variants <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7)

# extract useful columns
BeatAML_variants_sub <- BeatAML_variants %>%
  select("labId", "symbol", "t_vaf", "variant_class", "short_aa_change")

colnames(BeatAML_variants_sub)[4:5] <- c("variant_type", "amino_acid_change")

# remove duplicate calls for the same mutation but different VAFs
BeatAML_variants_sub <- BeatAML_variants_sub %>% 
  distinct(labId, symbol, t_vaf, amino_acid_change, .keep_all = T)


# subset to patients with required data types
# select pts with mutation data
pt_subset_1 <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$exomeSeq == "y")

#select patients with sequencing performed on BM aspirates
pt_subset_bm <- subset(pt_subset_1, pt_subset_1$specimenType == "Bone Marrow Aspirate")
pt_subset_bm <- pt_subset_bm %>% distinct(PatientId, .keep_all = TRUE)

# find unique pts with assays performed on peripheral blood that are not duplicates on the bone marrow assays
pt_subset_pb <- subset(pt_subset_1, pt_subset_1$specimenType == "Peripheral Blood")
pt_subset_pb <- pt_subset_pb %>% distinct(PatientId, .keep_all = TRUE)
# pt_subset_pb$PatientId <- as.character(pt_subset_pb$PatientId)

pt_subset_pb <- setDT(pt_subset_pb)[!(PatientId) %in% pt_subset_bm$PatientId]

pt_subset_2 <- rbind(pt_subset_bm, pt_subset_pb)

pt_subset_2 <- pt_subset_2 %>% distinct(PatientId, .keep_all = TRUE)
pt_subset_2 <- select(pt_subset_2, LabId, PatientId, isDenovo, isRelapse, isTransformed)
pt_subset_2$isDenovo <- as.character(pt_subset_2$isDenovo)
pt_subset_2$isRelapse <- as.character(pt_subset_2$isRelapse)
pt_subset_2$isTransformed <- as.character(pt_subset_2$isTransformed)
pt_subset_2 <- na.omit(pt_subset_2)

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

pt_subset_2 <- select(pt_subset_2, LabId, Subset)
colnames(pt_subset_2)[1] <- "labId"

# filter to variants present in the patient group and ensure no duplicate pts are represented
BeatAML_variants_sub <- setDT(BeatAML_variants_sub)[labId %chin% pt_subset_2$labId] 
BeatAML_variants_sub <- BeatAML_variants_sub %>%
  distinct(labId, symbol, .keep_all = TRUE)

BeatAML_variants_sub <- unique(BeatAML_variants_sub)
colnames(BeatAML_variants_sub)[3] <- "VAF"
combined <- BeatAML_variants_sub
combined$VAF <- as.numeric(combined$VAF) 

colnames(combined)[4] <- c("variant_type")

combined <- left_join(BeatAML_variants_sub, pt_subset_2, by = "labId")

# replace the labId with Patient ID in order to leverage the survival data
key <- select(BeatAML_sample_data_types, LabId, PatientId)
colnames(key)[1] <- "labId"

key <- setDT(key)[labId %chin% combined$labId] 

combined <- left_join(key, combined, by = "labId")

# remove labid column
combined$labId <- NULL
colnames(combined)[1] <- "labId"


### combine the BeatAML and TCGA mutation files together to find mutation frequencies and co-occurences
mut_table_final_TCGA_sub <- mut_table_final_TCGA[,c(1,2,5,3,6,9)]

# homogenize variant type annotations between the two datasets
for (i in 1:nrow(mut_table_final_TCGA_sub)) {
  if(mut_table_final_TCGA_sub$type[i] == "DEL"){
    mut_table_final_TCGA_sub$type[i] <- "deletion"
  }
  if(mut_table_final_TCGA_sub$type[i] == "INS"){
    mut_table_final_TCGA_sub$type[i] <- "insertion"
  }
  if(mut_table_final_TCGA_sub$type[i] == "SNP"){
    mut_table_final_TCGA_sub$type[i] <- "SNV"
  }
}


# make sure the column headers are the same for rowbind
cnames <- colnames(combined)
colnames(mut_table_final_TCGA_sub) <- cnames

# remove duplicate rows in both column prior to merging
mut_table_final_BeatAML_sub <- unique(BeatAML_variants_sub)
mut_table_final_TCGA_sub <- unique(mut_table_final_TCGA_sub)
colnames(mut_table_final_TCGA_sub)[3] <- "VAF"
mut_table_final_TCGA_sub$VAF <- as.numeric(mut_table_final_TCGA_sub$VAF)



combined <- rbind(combined, mut_table_final_TCGA_sub)
colnames(combined)[4] <- "variant_type"

#
#
#
#
#
#
#


#### Stanford: add samples to the other cohorts ####
# load list of file names and rename column header
files_list <- as.data.frame(list.files(path="~/Desktop/MetaAML_results/raw_data/Archive/", pattern="*.txt", full.names=F, recursive=FALSE))

colnames(files_list)[1] <- "file_path"

# create empty list to populate results into 
results_list <- list()
z <- 1

# loop through the data frame of individual file names and extract the individual mutations and VAFs
for(i in 1:nrow(files_list)){
  
  # i <- 7
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
    # print("here2")
    
    results_list[[z]] <- temp_dat_final
    z <- z + 1 
  } else {
    results_list[[z]] <- temp_dat1
    
    z <- z + 1 
  }
}


# bind the list of lists together to generate the final data frame
temp_final = as.data.frame(do.call(rbind, results_list))


clin_annotations <- read_excel("~/Desktop/MetaAML_results/raw_data/Archive/Copy of 160104_RM_AML_Sample_Info.xlsx")

# select columns of interest
labels <- select(clin_annotations, Sample, `Disease Status`)
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
    labels$`Disease Status`[i] <- "transformed"
  }
  if(labels$`Disease Status`[i] == "?"){
    labels$`Disease Status`[i] <- "other"
  }
  if(labels$`Disease Status`[i] == "Unknown"){
    labels$`Disease Status`[i] <- "other"
  }
}

colnames(labels)[2] <- "Subset"


clin_annotations <- select(clin_annotations, Sample, `Clonal Mutations (genotyped by Ryan)`)


flt3_itd <- clin_annotations[grepl("FLT3-ITD", clin_annotations$`Clonal Mutations (genotyped by Ryan)`), ]
colnames(flt3_itd)[2] <- "Mutation"
flt3_itd$Mutation <- "FLT3"
flt3_itd$variant_type <- "ITD"


temp_final_flt3 <- left_join(temp_final, flt3_itd, by = c("Sample" = "Sample", "Mutation" = "Mutation"))
# temp_final_flt3$VAF <- NULL
temp_final_flt3 <- unique(temp_final_flt3)

temp_final_flt3[is.na(temp_final_flt3)] <- 0

# label the mutation as FLT3-ITD or FLT3-TKD
for(i in 1:nrow(temp_final_flt3)){
  # print(i)
  if(temp_final_flt3$variant_type[i] == "ITD"){
    temp_final_flt3$SNV_or_Indel[i] <- "ITD"
  }
}

# remove unnecessary column  
temp_final_flt3$variant_type <- NULL

colnames(temp_final_flt3)[1:4] <- c("labId", "symbol", "VAF", "variant_type")

### seems like there are rare mutations in the stanford samples that are not recurrent in the other studdies. For now, filter to only mutations found in the other studdies
temp_final_flt3 <- setDT(temp_final_flt3)[symbol %chin% variants$gene_name] 

temp_final_flt3$VAF <- sub("%$", "", temp_final_flt3$VAF)
temp_final_flt3$VAF <- as.numeric(as.character(temp_final_flt3$VAF))

temp_final_flt3$VAF <- (temp_final_flt3$VAF/100)

# add the sample type column
temp_final_flt3 <- left_join(temp_final_flt3, labels, by = "labId")

# create the combined dataframe for BeatAML, TCGA, and Stanford
combined_2 <- rbind(combined, temp_final_flt3)


#
#
#
#
#
#
#



## Multistage: add the AML knowledgd bank mutation data ####
kb_data_muts <- read.table(text = x, header = T, stringsAsFactors = FALSE)

# use flt3 calls to annotation mutations according to type
kb_data_muts_sub <- select(kb_data_muts, SAMPLE_NAME, VAF, VARIANT_TYPE, GENE, AA_CHANGE)

clinicalData2 <- read.table("~/Desktop/MetaAML_results/raw_data/AML_knowledge_bank_data_clinical.txt", stringsAsFactors = F)

# save data for survival analysis
Multistage_survival <- clinicalData2

clin_annotations2 <- select(clinicalData2, PDID, TypeAML, FLT3_ITD, FLT3_TKD)
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

# # kb_data_flt3 <- kb_data_flt3[,1:4]
# 
# names <- colnames(kb_data_flt3)
# 
# joined_muts_clinical <- left_join(kb_data_muts_sub, kb_data_flt3,  by = "SAMPLE_NAME")
# 
# joined_muts_clinical <- joined_muts_clinical[,c(1,5:11)]
# colnames(joined_muts_clinical) <- c(names)
joined_muts_clinical <- kb_data_flt3

# rename the varient type for homogeneity with the other data sets
for(i in 1:nrow(joined_muts_clinical)){
  
  if(joined_muts_clinical$VARIANT_TYPE[i] == "D" & joined_muts_clinical$FLT3_ITD[i] == 1 & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "I" & joined_muts_clinical$FLT3_ITD[i] == 1 & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "ID" & joined_muts_clinical$FLT3_ITD[i] == 1 & joined_muts_clinical$GENE[i] == "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "D" & joined_muts_clinical$GENE[i] != "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "deletion"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "I" & joined_muts_clinical$GENE[i] != "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "insertion"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "ID" & joined_muts_clinical$GENE[i] != "FLT3"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "Sub"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "SNV"
  }
  if(joined_muts_clinical$VARIANT_TYPE[i] == "PTD"){
    joined_muts_clinical$VARIANT_TYPE[i] <- "ITD"
  }
}

joined_muts_clinical[,7:8] <- NULL
joined_muts_clinical <- unique(joined_muts_clinical)

joined_muts_clinical <- select(joined_muts_clinical, SAMPLE_NAME, GENE, VARIANT_TYPE, VAF, AA_CHANGE, TypeAML)
colnames(joined_muts_clinical)[1:6] <- c("labId", "symbol", "variant_type", "VAF", "amino_acid_change", "Subset")

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


# normalize the VAF to decimal format
joined_muts_clinical$VAF <- as.numeric(joined_muts_clinical$VAF, na.rm = T)
joined_muts_clinical$VAF <- (joined_muts_clinical$VAF/100)


## make final combined dataframe
final_data_matrix_sub <- rbind(combined_2, joined_muts_clinical)
# make sure to make the VAF in the appropriate range
final_data_matrix_sub$VAF <- round(final_data_matrix_sub$VAF, 2)



## for some reason, SRSF2 used to be SFRS2 and was not changed in any of the datasets except BeatAML. Find the old annotations and change them here:
for(i in 1:nrow(final_data_matrix_sub)){
  if(final_data_matrix_sub$symbol[i] == "SFRS2"){
    final_data_matrix_sub$symbol[i] <- "SRSF2"
  }
}

colnames(final_data_matrix_sub)[1] <- "Sample"



#### Ley et al 2012 ####
# select necessary columns in the mutation data
ley_mut_2_sub <- select(ley_mut, AML, 'Annotated Gene', 'Mutation Type', 'Variant Allele Frequency', 'Mutation Effect')
ley_mut_2_sub$Subset <- "de_novo"
colnames(ley_mut_2_sub) <- c("labId", "symbol", "variant_type", "VAF", "amino_acid_change", "Sample")

#
#
#
#
#
#


## Now include the survival data to the table

TCGA_survival2 <- select(TCGA_survival, PATIENT_ID, OS_STATUS, OS_MONTHS, SEX, CYTOGENETIC_CODE_OTHER, RISK_MOLECULAR)
# need to convert tcga data from months to days
TCGA_survival2$OS_MONTHS <- (TCGA_survival2$OS_MONTHS*30)

BeatAML_survival2 <- select(BeatAML_sample_data_types, PatientId, vitalStatus, overallSurvival, consensus_sex, Karyotype, ELN2017)
Stanford_survival2 <- select(Stanford_survival, `Patient SU Number`, `Overall Survival (Event = Death)`, `Duration from diagnosis to D1`, Gender, Cytogenetics, 15)
Stanford_survival2 <- Stanford_survival2[1:133,1:6]

Multistage_survival2 <- select(Multistage_survival, PDID, Status, OS, gender, complex, M_Risk)

# make all the headers the same
# mutation_matrix <- read.csv("~/Desktop/Majeti_Lab/Data/Combined_mutation_occurence/De_novo/MetaAML_mutation_matrix_De_novo_15_pts.csv", header = T)
# mutation_matrix$X <- NULL

colnames(TCGA_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Cytogenetics", "Risk")
colnames(BeatAML_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Cytogenetics", "Risk")
colnames(Stanford_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Cytogenetics", "Risk")
colnames(Multistage_survival2) <- c("Sample", "Censor", "Time_to_OS", "Sex", "Cytogenetics", "Risk")

mutation_matrix_survival <- rbind(TCGA_survival2, BeatAML_survival2, Stanford_survival2, Multistage_survival2)
mutation_matrix_survival <- na.omit(mutation_matrix_survival)
# print("here3")


for(i in 1:nrow(mutation_matrix_survival)){
  if(mutation_matrix_survival$Censor[i] == "LIVING" | mutation_matrix_survival$Censor[i] == "Alive"){
    mutation_matrix_survival$Censor[i] <- 1
  }
  if(mutation_matrix_survival$Censor[i] == "DECEASED" | mutation_matrix_survival$Censor[i] == "Dead"){
    mutation_matrix_survival$Censor[i] <- 0
  }
}

# remove cases where the time to OS is negative
mutation_matrix_survival <- mutation_matrix_survival[mutation_matrix_survival$Time_to_OS >= 0, ]

final_data_matrix <- left_join(final_data_matrix_sub, mutation_matrix_survival,  by = "Sample")

# homogenize the different terms for each dataset
for(i in 1:nrow(final_data_matrix)){
  if(!is.na(final_data_matrix$Sex[i])){
    if(final_data_matrix$Sex[i] == "M" | final_data_matrix$Sex[i] == "1" | final_data_matrix$Sex[i] == "Male"){
      final_data_matrix$Sex[i] <- "Male"
    }
    if(final_data_matrix$Sex[i] == "F" | final_data_matrix$Sex[i] == "2" | final_data_matrix$Sex[i] == "Female"){
      final_data_matrix$Sex[i] <- "Female"
    } 
  }
  if(!is.na(final_data_matrix$Risk[i])){
    if(final_data_matrix$Risk[i] == 0 | final_data_matrix$Risk[i] == "Good"){
      final_data_matrix$Risk[i] <- "Favorable"
    } 
    if(final_data_matrix$Risk[i] == 1 | final_data_matrix$Risk[i] == "Inter-1" | final_data_matrix$Risk[i] == "Inter-2"){
      final_data_matrix$Risk[i] <- "Intermediate"
    } 
    if(final_data_matrix$Risk[i] == 2 | final_data_matrix$Risk[i] == "Poor" | final_data_matrix$Risk[i] == "Inter-2"){
      final_data_matrix$Risk[i] <- "Adverse"
    } 
    if(final_data_matrix$Risk[i] == "NA" | final_data_matrix$Risk[i] == "Unknown" | final_data_matrix$Risk[i] == "N.D."){
      final_data_matrix$Risk[i] <- NA
    } 
  }
}


# find the different annotation for cytogenetics
# unique_cyto <- as.data.frame(unique(final_data_matrix$Cytogenetoics))
final_data_matrix <- final_data_matrix[!duplicated(final_data_matrix[1:5]),]

save(final_data_matrix,  file = "~/Desktop/MetaAML_results/final_data_matrix.RData")

# subset for AMy
MetaAML_RUNX1 <- subset(final_data_matrix, final_data_matrix$symbol == "RUNX1")
write.csv(MetaAML_RUNX1,  file = "~/Desktop/Majeti_Lab/Data/MetaAML/RUNX1_data.csv", row.names=FALSE)
