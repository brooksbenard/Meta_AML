# MetaAML_oncoprint
#
# Brooks Benard
# bbenard@stanford.edu
# 10/09/2019
#

### This script summarizes the major genomics and  MetaAML to visualize the MetaAML cohort

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('ComplexHeatmap')) install.packages('ComplexHeatmap'); library('ComplexHeatmap')

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

# final_data_matrix_3 <- filter(final_data_matrix_3, VAF >= 0.01 & final_data_matrix_3$Subset == "de_novo")

final_data_matrix_3 <- subset(final_data_matrix_3, final_data_matrix_3$mut_freq_gene >= 100)

c = as.numeric(n_distinct(final_data_matrix_3$Sample))
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

temp_dat <- as.matrix(temp_dat)

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

 plot = oncoPrint(temp_dat, 
                  col = col,
                  # top_annotation_height = 2,
                  row_names_side= "left",
                  # show_pct = T,
                  bottom_annotation = ha,
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

 print(plot)
 
 png(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_oncoprint_100.png", res = 300, width = 15, height = 6.5, units = "in")
 
 plot
 print(plot)
 dev.off()
 
 
  # ggsave(plot = plot, filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_oncoprint.png", dpi = 300, width = 10, height = 5, units = "in")
  

