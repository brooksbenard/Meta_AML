# load required packages
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
if (!require('vegan')) install.packages('vegan'); library('vegan')


# download data from BeatAML website
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/Majeti_Lab/Data/BeatAML/41586_2018_623_MOESM3_ESM.xlsx")


BeatAML_sample_data_types <- read_excel("~/Desktop/Majeti_Lab/Data/BeatAML/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)

# subset to patients with required data types
pt_subset_1 <- subset(BeatAML_sample_data_types, BeatAML_sample_data_types$exomeSeq == "y" & BeatAML_sample_data_types$totalDrug == "y")


# select patient subset
  pt_subset_2 <- subset(pt_subset_1, pt_subset_1$isDenovo == "TRUE" & pt_subset_1$isRelapse == "FALSE")
  group <- "De_novo"
  label_1 <- "De novo"

# } else if(subset == "Relapse"){
#   pt_subset_2 <- subset(pt_subset_1, pt_subset_1$isRelapse == "TRUE")
#   group <- "Relapse"
#   label_1 <- "Relapse"
#   
# } else if(subset == "Transformed"){
#   pt_subset_2 <- subset(pt_subset_1, pt_subset_1$isTransformed == "TRUE" & pt_subset_1$isRelapse == "FALSE")
#   group <- "Transformed"
#   label_1 <- "Secondary"
# } else if(subset == "All"){
#   pt_subset_1 <- pt_subset_1[!duplicated(pt_subset_1$LabId), ]
#   
#   # find bm samples
#   pt_subset_2 <- subset(pt_subset_1, pt_subset_1$specimenType == "Bone Marrow Aspirate")
#   
#   # find peripheral blood samples and filter to only those without bm samples
#   pt_subset_3 <- subset(pt_subset_1, pt_subset_1$specimenType == "Peripheral Blood")
#   pt_subset_3 <- setDT(pt_subset_3)[! PatientId %chin% pt_subset_2$PatientId] 
#   
#   pt_subset_2 <- rbind(pt_subset_2, pt_subset_3)
#   
#   group <- "All"
#   label_1 <- "All"
#   
# }


# read in BeatAML variants for analysis data
BeatAML_variants <- read_excel("~/Desktop/Majeti_Lab/Data/BeatAML/41586_2018_623_MOESM3_ESM.xlsx", sheet = 7)

BeatAML_variants <- BeatAML_variants %>%
  select("labId", "symbol", "t_vaf", "short_aa_change", "chrom")

# deliniate FLT3 from FLT3-ITD
for(i in 1:nrow(BeatAML_variants)){
  print(i)
  if(BeatAML_variants[i,2] == "FLT3" & BeatAML_variants[i,4] == "ITD"){
    BeatAML_variants[i,2] <- "FLT3_ITD"
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


# add concensus sex column
gender <- pt_subset_2[,1:3]
colnames(gender)[1] <- "labId"
BeatAML_variants_sub <- dplyr::left_join(BeatAML_variants_sub, gender, by = "labId")

# for(i in 1:nrow(BeatAML_variants_sub)){
#   print(i)
#   if(BeatAML_variants_sub$consensus_sex[i] == "Male" && BeatAML_variants_sub$chrom[i] == "X"){
#     BeatAML_variants_sub$t_vaf[i] <- (BeatAML_variants_sub$t_vaf[i])/2
#   }
# }

mut_table <- aggregate(data.frame(count = BeatAML_variants_sub), list(value = BeatAML_variants_sub$symbol), length)


# read in drug response data
drug_data <- read_excel("~/Desktop/Majeti_Lab/Data/BeatAML/41586_2018_623_MOESM3_ESM.xlsx", sheet = 10)

## pair patient samples with drug response data
# filter to only patients with both mutation and drug screening data
colnames(BeatAML_variants_sub)[1] <- "lab_id"
drug_data <- setDT(drug_data)[lab_id %chin% BeatAML_variants_sub$lab_id] 

drug_mut <- dplyr::right_join(drug_data, BeatAML_variants_sub, by = "lab_id")


## format drug list
# they have one drug that is not formated like all the others so I manually reformat it before the loop
drug_mut[] <- lapply(drug_mut, gsub, pattern = "17-AAG (Tanespimycin)", replacement = "Tanespimycin", fixed = TRUE)

# loop through drugs and create a matrix of all drug-mutation-vaf combinations
drug_list <- list()
for(i in 1:nrow(drug_mut)){
  print(i)
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



# calculate shannon diversity index ####
samples=as.data.frame(unique(drug_mut$lab_id))
auc_sd_list = list()

z = 1

for(i in 1:nrow(samples)){
  print(i)
  sub = subset(drug_mut, drug_mut$lab_id == samples[i,1])
  sub = distinct(sub,lab_id, t_vaf)
  
  if(nrow(sub) > 1){
    
    sd=diversity(x = sub$t_vaf, index = "shannon")
    
    auc_sd <- data.frame(matrix(NA, nrow = 1, ncol = 2))
    names(auc_sd) <- c("lab_id", "shannon_diversity_index")
    
    auc_sd[1,1] <- sub[1,1]
    auc_sd[1,2] <- sd
    
    # Add each list in the loop to a list of lists
    z = z + 1
    auc_sd_list[[z]] <- auc_sd  
  }
}

beatAML_auc_shannon_diversity <- do.call(rbind, auc_sd_list)


# append to drug results dataframe