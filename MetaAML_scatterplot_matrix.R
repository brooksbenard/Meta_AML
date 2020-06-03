# Brooks Benard
# bbenard@stanford.edu
# 03/23/2020
#

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('patchwork')) install.packages('patchwork'); library('patchwork')
if (!require('rlist')) install.packages('rlist'); library('rlist')

load("~/Desktop/MetaAML_results/final_data_matrix.RData")

scatterplot_data = subset(final_data_matrix, final_data_matrix$mut_freq_gene >= 150 & final_data_matrix$Gene != "MLL" & final_data_matrix$Subset == "de_novo")

# make sure that the FLT3 symbols are annotated properly
for(i in 1:nrow(scatterplot_data)){
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "ITD"){
    scatterplot_data$Gene[i] <- "FLT3-ITD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "SNV"){
    scatterplot_data$Gene[i] <- "FLT3-TKD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "Deletion"){
    scatterplot_data$Gene[i] <- "FLT3-TKD"
  }
  if(scatterplot_data$Gene[i] == "FLT3" & scatterplot_data$variant_type[i] == "INDEL"){
    scatterplot_data$Gene[i] <- "FLT3-ITD"
  }
}


# select useful columns
scatterplot_data = select(scatterplot_data, Sample, Gene, VAF)

# remove FLT3-ITD rows because they have no VAF data
scatterplot_data = subset(scatterplot_data, Gene != "FLT3-ITD")

# function to make a scatterplot for each pairwise genes entered
vaf_scatterplot_function <- function(gene_1, gene_2){
  
  # select mutations of interest
  sub_gene_1 <- subset(scatterplot_data, scatterplot_data$Gene == gene_1)
  sub_gene_2 <- subset(scatterplot_data, scatterplot_data$Gene == gene_2)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  if(nrow(gene_1_and_2) > 4){
  
      gene_1_and_2=as.data.frame(gene_1_and_2) %>% distinct(Sample, VAF.x, VAF.y, .keep_all = F)
    
    gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF.x),]
    gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
    
    # add point color columns for visualizing clonal/subclonal trends
    gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF.x - gene_1_and_2$VAF.y)
    
    gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
    
    if(nrow(gene_1_and_2) > 4){
      
      # define point color
      gene_1_and_2$Clonality <- NA
      gene_1_and_2$Clonality_color <- NA
      
      for(i in 1:nrow(gene_1_and_2)){
        if(gene_1_and_2$vaf_ratio[i] <= .05 & gene_1_and_2$vaf_ratio[i] >= -.05){
          gene_1_and_2$Clonality[i] <- "Ambiguous"
          gene_1_and_2$Clonality_color[i] = "Ambiguous"
        }
        if(gene_1_and_2$vaf_ratio[i] > .05){
          gene_1_and_2$Clonality[i] <- paste(gene_1, "first", sep = " ")
          gene_1_and_2$Clonality_color[i] = "gene1_first"
        }
        if(gene_1_and_2$vaf_ratio[i] < -.05){
          gene_1_and_2$Clonality[i] <- paste(gene_2, "first", sep = " ")
          gene_1_and_2$Clonality_color[i] = "gene2_first"
        }
      }
      
      gene_1_and_2$Clonality_color <- as.factor(gene_1_and_2$Clonality_color)
      
      
      # if(nrow(gene_1_and_2) > 4){
        # make the scatterplot
        scatter_plot = ggplot(gene_1_and_2, aes(x = gene_1_and_2$VAF.y, y = gene_1_and_2$VAF.x)) +
          geom_point(aes(color = Clonality_color), size =3, alpha = 0.75) +
          scale_x_continuous(
            limits = c(0,1),
            labels = scales::number_format(accuracy = 0.1)) +
          scale_y_continuous(
            limits = c(0,1),
            labels = scales::number_format(accuracy = 0.1)) +
          xlab(paste(gene_2, " VAF", sep ="")) +
          ylab(paste(gene_1, " VAF", sep ="")) +  
          scale_color_manual(values = c("Ambiguous" = "#374E55FF", "gene2_first" = "#8c510a", "gene1_first" = "#01665e")) +
          geom_abline(intercept = .05, slope = (1), color="#969696",
                      linetype="dashed", size=.5)+
          geom_abline(intercept = -.05, slope = (1), color="#969696",
                      linetype="dashed", size=.5)+
          geom_point(shape = 1, size =  3,colour = "black") +
          theme(legend.position = "none", 
                legend.title = element_blank()
                ) 
      # }    
    }
  }
}

# find all pairwise gene comparisons
scatterplot_genes=as.data.frame(t(combn(as.vector(unique(scatterplot_data$Gene)),2)))

scatterplot_genes <- scatterplot_genes[order(scatterplot_genes$V1),]

# make list to populate plots into 
sc_plots = list()

for(j in 1:nrow(scatterplot_genes)){
  print(j)
  g1 = as.character(scatterplot_genes[j,1])
    g2 = as.character(scatterplot_genes[j,2])
    
  sc_plots[[j]] = vaf_scatterplot_function(gene_1 = g1, gene_2 = g2)
}

sc_plots = list.clean(sc_plots)

p = cowplot::plot_grid(plotlist = sc_plots, 
                       ncol = 7, align = "hv",   
                       rel_heights = c(1,1),
                       rel_widths = c(1,1))

ggsave(filename ="~/Desktop/MetaAML_results/Data/Figures/vaf_scatterplots/MetaAML_scatterplot_matrix.png", dpi = 150, width = 17, height = 22, units = "in")

p = cowplot::plot_grid(plotlist = sc_plots, 
                       ncol = 13, align = "hv",   
                       rel_heights = c(1,1),
                       rel_widths = c(1,1))

ggsave(filename ="~/Desktop/MetaAML_results/Data/Figures/vaf_scatterplots/MetaAML_scatterplot_matrix_ppt.png", dpi = 150, width = 25, height = 13, units = "in")

rm(sc_plots)
