### MetaAML_data_upset_plot
## Brooks Benard
## bbenard@stanford.edu
## 02092020
# 
# This script plots the overlap of different data types for the meta analysis cohort
#
# load required packages
if (!require('UpSetR')) install.packages('UpSetR'); library('UpSetR')

load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_2 = final_data_matrix

# download annotation data from BeatAML website
# download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/Majeti_Lab/Data/BeatAML/41586_2018_623_MOESM3_ESM.xlsx")

# annotate RNA-seq and drug screened samples
BeatAML_sample_data_types <- read_excel("~/Desktop/MetaAML_results/raw_data/41586_2018_623_MOESM3_ESM.xlsx", sheet = 5)
BeatAML_sample_data_types = select(BeatAML_sample_data_types, LabId, PatientId, totalDrug, rnaSeq)
BeatAML_sample_data_types$rna_seq = ifelse(BeatAML_sample_data_types$rnaSeq == "y", 1,0)
BeatAML_sample_data_types$drug_screening = ifelse(BeatAML_sample_data_types$totalDrug == "y", 1,0)
BeatAML_sample_data_types = select(BeatAML_sample_data_types, LabId, rna_seq, drug_screening)
colnames(BeatAML_sample_data_types)[1] = "Sample"


# TCGA expression data annotations
# download.file("https://api.gdc.cancer.gov/data/57b0d7aa-5e00-4deb-86ac-9f22f819df9c", destfile = "~/Desktop/Majeti_Lab/Data/CD123/laml_tcga_pub.tar.gz")
# untar("~/Desktop/Majeti_Lab/Data/CD123/laml_tcga_pub.tar.gz")
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
subset_data = select(final_data_matrix_2, Sample, Gene, VAF, Subset, Cohort, Time_to_OS)

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
pdf(file = "~/Desktop/MetaAML_results/Data/Figures/cohort_upset_plot.pdf", width = 5.5, height = 4)
p
dev.off()
