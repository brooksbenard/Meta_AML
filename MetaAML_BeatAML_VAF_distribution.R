### This script plots the VAF disctibution for recurrent AML mutations in the BeatAML cohort
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
# download data from BeatAML website
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx", destfile = "~/Desktop/Majeti_Lab/Data/BeatAML/41586_2018_623_MOESM3_ESM.xlsx")

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
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF") +
  geom_jitter(aes(fill = variant_class), color = "black", shape = 21, position=position_jitter(0.2), size = 2) +
  scale_fill_manual(values = c("#6A6599FF", "#DF8F44FF","#79AF97FF","#B24745FF")) +
  # geom_jitter(shape=21, position=position_jitter(0.2)) +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label = "VAF") +
  xlab(label = NULL) +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))

ggpar(p, legend.title = "Variant class")

ggsave(filename = "~/Desktop/MetaAML_results/Data/drug_vaf_correlation/vaf_distribution_de_vovo.pdf", dpi = 300, width = 9, height = 4.5, units = "in")



# MetaAML VAf distribution summary ####
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

vaf_sub$threshold = ifelse(vaf_sub$VAF_male_x >= .3 , "Clonal","Subclonal")

vaf_sub = subset(vaf_sub, vaf_sub$threshold == "Clonal" | vaf_sub$threshold == "Subclonal")

vaf_sub$VAF_male_x=as.numeric(vaf_sub$VAF_male_x)
vaf_sub$Gene = as.character(vaf_sub$Gene)

vaf_sub$Gene <- with(vaf_sub, reorder(Gene, -VAF_male_x, median))


p = ggplot(vaf_sub, aes(x=Gene, y=VAF_male_x)) + 
  geom_boxplot(notch=F, outlier.colour = "white", color = "#374E55FF", fill = "lightgrey") +
  geom_jitter(aes(fill = threshold), color = "black", shape = 21, position=position_jitter(0.2), size = 1.5) +
  scale_fill_manual(values = c("#cb181d", "#3690c0")) +
  # scale_fill_manual(values = c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF", SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")) +
  geom_hline(yintercept = .3, color = "#b2182b", linetype="dashed") +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label = "VAF") +
  xlab(label = NULL) +
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))

ggpar(p, legend.title = "")
ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/MetaAML_vaf_distribution.pdf", dpi = 300, width = 10, height = 4, units = "in")
