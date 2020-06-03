# MATH (Mutant-Allele Tymor Heterogeneity) score calculation ####
# calculate an additional number for tumor heterogeneity for each patient
# MATH = 100 * MAD/median

pts=as.data.frame(unique(final_data_matrix_2$Sample))

colnames(pts)[1]="ID"
pts$ID=as.character(pts$ID)

math_score_list_2=list()
z=1
for(i in 1:nrow(pts)){
  print(i)
  pt=as.character(pts[i,1])
  sub=subset(final_data_matrix_2, final_data_matrix_2$Sample==pt)
  if(nrow(sub > 1)){
    # MATH = 100 * MAD/median
    abs.med.dev = abs(sub$VAF_male_x - median(sub$VAF_male_x)) #absolute deviation from median vaf
    pat.mad = median(abs.med.dev) * 100
    pat.math = pat.mad * 1.4826 /median(sub$VAF_male_x)
    
    # mad_ = 1.4826*(mad(sub$VAF_male_x))    
    # median_ = median(sub$VAF_male_x)
    # math_ = 100*(mad_/median_)
    
    math_score_list <- data.frame(matrix(NA, nrow = 1, ncol = 2))
    names(math_score_list) <- c("Sample", "MATH_score")
    
    math_score_list[1,1] <- pt
    math_score_list[1,2] <- pat.math
    
    # Add each list in the loop to a list of lists
    z <- z + 1
    math_score_list_2[[z]] <- math_score_list 
  }
}
math_score_list_final = as.data.frame(do.call(rbind, math_score_list_2))

# append the shannon diversity index to the dataframe
final_data_matrix_2=left_join(final_data_matrix_2, math_score_list_final, by = "Sample")

# MATH by number of mutations ####
final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")

# for visualization purposes, collapse the patients with more than 7 mutations for plotting the shannon diversity index
final_data_matrix_2_sub$mut_freq_pt_collapsed = NA
for(i in 1:nrow(final_data_matrix_2_sub)){
  if(final_data_matrix_2_sub$mut_freq_pt[i] > 7){
    final_data_matrix_2_sub$mut_freq_pt_collapsed[i] = "8+"
  }
  if(final_data_matrix_2_sub$mut_freq_pt[i] <= 7){
    final_data_matrix_2_sub$mut_freq_pt_collapsed[i] = final_data_matrix_2_sub$mut_freq_pt[i]
  }
}

final_data_matrix_2_sub = subset(final_data_matrix_2_sub, final_data_matrix_2_sub$mut_freq_pt > 2)

final_data_matrix_2_sub <- na.omit(distinct(final_data_matrix_2_sub, Sample, mut_freq_pt_collapsed, MATH_score))

# # pal=pal_uchicago()(8)
ggplot(final_data_matrix_2_sub, aes(x=mut_freq_pt_collapsed, y=MATH_score, fill = mut_freq_pt_collapsed)) + 
  geom_boxplot(notch=F, outlier.colour = "white") +
  # stat_boxplot(geom ='errorbar') +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  # scale_x_discrete(limits=c("2", "3", "4", "5", "6", "7" ,"8+")) +
  # stat_compare_means() +
  scale_fill_uchicago() +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label= "MATH score") +
  xlab(label = "Number of unique mutations per patient") +
  theme(legend.position="none") 

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/MATH_score_by_mutation_burden_de_novo.pdf", dpi = 300, width = 6, height = 4, units = "in")
