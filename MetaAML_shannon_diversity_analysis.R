# MetaAML_shannon_diversity_analysis
#
# Brooks Benard
# bbenard@stanford.edu
# 03/24/2020
#
#### load required packages and data ####
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis'); library('survivalAnalysis')
if (!require('vegan')) install.packages('vegan'); library('vegan')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')

# load data
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# Shannon diversity calcluation ####
# calculate a diversity index based on the the vafs for each patient
pts=as.data.frame(unique(final_data_matrix_2$Sample))
colnames(pts)[1]="ID"
pts$ID=as.character(pts$ID)

sub = select(final_data_matrix_2, Sample, Gene, VAF_male_x, clonality_bin)

sd_list2=list()
z=1
for(i in 1:nrow(pts)){
  # print(i)
  pt=as.character(pts[i,1])
  sub_2 = na.omit(subset(sub, sub$Sample==pt))
  
  # ensure there is enough data to perform the diverity analysis
  if(n_distinct(sub_2$clonality_bin) > 1){
    
    # estimate cancer cell fractions from the average vaf in each "clonal" bin
    ccf <- aggregate(sub_2$VAF_male_x, by = list(sub_2$clonality_bin), mean)
    
    sd <- diversity(x = ccf$x, index = "shannon")
    
    sd_list <- data.frame(matrix(NA, nrow = 1, ncol = 2))
    names(sd_list) <- c("Sample", "Shannon_diversity")
    
    sd_list[1,1] <- pt
    sd_list[1,2] <- sd
    
    # Add each list in the loop to a list of lists
    z <- z + 1
    sd_list2[[z]] <- sd_list 
  }
}
sd_list_final = as.data.frame(do.call(rbind, sd_list2))

# append the shannon diversity index to the dataframe
final_data_matrix_2=left_join(final_data_matrix_2, sd_list_final, by = "Sample")


save(final_data_matrix_2,  file = "~/Desktop/MetaAML_results/final_data_matrix_2.RData")
rm(list=setdiff(ls(), c("final_data_matrix", "final_data_matrix_2")))






# Shannon Diversity by number of mutations ####
# final_data_matrix_2_sub <- subset(final_data_matrix_2, !Cohort %in% c("Au","Wang","Garg","Huet"))
# final_data_matrix_2_sub <- subset(final_data_matrix_2_sub, Gene!="MLL")

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

# final_data_matrix_2_sub = subset(final_data_matrix_2_sub, final_data_matrix_2_sub$mut_freq_pt != 1)

final_data_matrix_2_sub <- na.omit(distinct(final_data_matrix_2_sub, Sample, mut_freq_pt_collapsed, Shannon_diversity))

# # pal=pal_uchicago()(8)
ggplot(final_data_matrix_2_sub, aes(x=mut_freq_pt_collapsed, y=Shannon_diversity, fill = mut_freq_pt_collapsed)) + 
  geom_boxplot(notch=F, outlier.colour = "white") +
  stat_boxplot(geom ='errorbar') +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  scale_x_discrete(limits=c("2", "3", "4", "5", "6", "7" ,"8+")) +
  # stat_compare_means() +
  scale_fill_uchicago() +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label= "Shannon Diversity Index") +
  xlab(label = "Number of unique mutations per patient") +
  theme(legend.position="none") 

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/Shannon_diversity_by_mutation_burden_de_novo.pdf", dpi = 300, width = 6, height = 4, units = "in")


# plot SD by age
final_data_matrix_2_sub <- final_data_matrix_2
final_data_matrix_2_sub <- subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
final_data_matrix_2_sub$age_bin <- NA

for(i in 1:nrow(final_data_matrix_2_sub)){
  if(!is.na(final_data_matrix_2_sub$Age[i])){
    if(final_data_matrix_2_sub$Age[i] <= 30){
      final_data_matrix_2_sub$age_bin[i] <- "< 30"
    }
    if(final_data_matrix_2_sub$Age[i] > 30 & final_data_matrix_2_sub$Age[i] <= 40){
      final_data_matrix_2_sub$age_bin[i] <- "30-40"
    }
    if(final_data_matrix_2_sub$Age[i] > 40 & final_data_matrix_2_sub$Age[i] <= 50){
      final_data_matrix_2_sub$age_bin[i] <- "40-50"
    }
    if(final_data_matrix_2_sub$Age[i] > 50 & final_data_matrix_2_sub$Age[i] <= 60){
      final_data_matrix_2_sub$age_bin[i] <- "50-60"
    }
    if(final_data_matrix_2_sub$Age[i] > 60 & final_data_matrix_2_sub$Age[i] <= 70){
      final_data_matrix_2_sub$age_bin[i] <- "60-70"
    }
    if(final_data_matrix_2_sub$Age[i] > 70){
      final_data_matrix_2_sub$age_bin[i] <- "70+"
    }
  }   
}

final_data_matrix_2_sub = subset(final_data_matrix_2_sub, final_data_matrix_2_sub$mut_freq_pt != 1)

final_data_matrix_2_sub <- na.omit(distinct(final_data_matrix_2_sub, Sample, age_bin, Shannon_diversity))


ggplot(final_data_matrix_2_sub, aes(x=age_bin, y=Shannon_diversity, fill = age_bin)) + 
  geom_boxplot(notch=F, outlier.colour = "white") +
  # stat_boxplot(geom ='errorbar') +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  # scale_x_discrete(limits=c("2", "3", "4", "5", "6", "7" ,"8+")) +
  stat_compare_means() +
  scale_fill_jco() +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label= "Shannon Diversity Index") +
  xlab(label = "Age") +
  theme(legend.position="none") 

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/Shannon_diversity_by_age_de_novo.pdf", dpi = 300, width = 6, height = 4, units = "in")



# shannon diversity by number of clones
final_data_matrix_2_sub = subset(final_data_matrix_2, final_data_matrix_2$Subset == "de_novo")
final_data_matrix_2_sub <- na.omit(distinct(final_data_matrix_2_sub, Sample, mut_freq_pt_collapsed, Shannon_diversity, n_clones))
final_data_matrix_2_sub$n_clones = as.factor(final_data_matrix_2_sub$n_clones)
# # pal=pal_uchicago()(8)
ggplot(final_data_matrix_2_sub, aes(x=n_clones, y=Shannon_diversity, fill = n_clones)) +
  geom_boxplot(notch=F, outlier.colour = "white") +
  stat_boxplot(geom ='errorbar') +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  # scale_x_discrete(limits=c("2", "3", "4")) +
  # stat_compare_means() +
  scale_fill_uchicago() +
  theme_cowplot(font_size = 15) +
  labs(title = NULL) +
  ylab(label= "Shannon Diversity Index") +
  xlab(label = "Number of clones") +
  theme(legend.position="none")
# 
ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/Shannon_diversity_by_mutation_burden_de_novo.pdf", dpi = 300, width = 6, height = 4, units = "in")


# survival by shannon diversity index ####
sub_shannon=data.frame(unique(select(final_data_matrix_2, Sample, Time_to_OS, Censor, Shannon_diversity, Subset, mut_freq_pt)))

sub_shannon=subset(sub_shannon, sub_shannon$Subset == "de_novo" & sub_shannon$mut_freq_pt > 1)

sub_shannon$mut_freq_pt_collapsed = NA
for(i in 1:nrow(sub_shannon)){
  if(sub_shannon$mut_freq_pt[i] > 7){
    sub_shannon$mut_freq_pt_collapsed[i] = "8+"
  }
  if(sub_shannon$mut_freq_pt[i] <= 7){
    sub_shannon$mut_freq_pt_collapsed[i] = sub_shannon$mut_freq_pt[i]
  }
}

sub_shannon$Shannon_diversity = as.numeric(sub_shannon$Shannon_diversity)

# find the median shannon diversity for each of the mutation bins
sd_two_mean=subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == 2)
sd_two_mean = median(sd_two_mean$Shannon_diversity, na.rm = T)
sd_three_mean = as.data.frame(subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == 3))
sd_three_mean = median(sd_three_mean$Shannon_diversity, na.rm = T)
sd_four_mean=subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == 4)
sd_four_mean = median(sd_four_mean$Shannon_diversity, na.rm = T)
sd_five_mean=subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == 5)
sd_five_mean = median(sd_five_mean$Shannon_diversity, na.rm = T)
sd_six_mean=subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == 6)
sd_six_mean = median(sd_six_mean$Shannon_diversity, na.rm = T)
sd_seven_mean=subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == 7)
sd_seven_mean = median(sd_seven_mean$Shannon_diversity, na.rm = T)
sd_more_mean=subset(sub_shannon, sub_shannon$mut_freq_pt_collapsed == "8+")
sd_more_mean = median(sd_more_mean$Shannon_diversity, na.rm = T)

# loop through the shannon diversity for each patient and assign if the patient is above or below the median sd for that mutaiton burden
sub_shannon$mean_sd = NA
sub_shannon$high_sd_or_low_sd = NA

# assign median sd for each mutation burden
for(i in 1:nrow(sub_shannon)){
  if(sub_shannon$mut_freq_pt[i] == 2){
    sub_shannon$mean_sd[i] = sd_two_mean
  }
  if(sub_shannon$mut_freq_pt[i] == 3){
    sub_shannon$mean_sd[i] = sd_three_mean
  }
  if(sub_shannon$mut_freq_pt[i] == 4){
    sub_shannon$mean_sd[i] = sd_four_mean
  }
  if(sub_shannon$mut_freq_pt[i] == 5){
    sub_shannon$mean_sd[i] = sd_five_mean
  }
  if(sub_shannon$mut_freq_pt[i] == 6){
    sub_shannon$mean_sd[i] = sd_six_mean
  }
  if(sub_shannon$mut_freq_pt[i] == 7){
    sub_shannon$mean_sd[i] = sd_seven_mean
  }
  if(sub_shannon$mut_freq_pt[i] > 7){
    sub_shannon$mean_sd[i] = sd_more_mean
  }
}

# assign if the patient is above or below the mean sd for theit mutation count
sub_shannon$high_sd_or_low_sd = ifelse(sub_shannon$Shannon_diversity > sub_shannon$mean_sd, 1,0)
sub_shannon$high_sd_or_low_sd = as.numeric(sub_shannon$high_sd_or_low_sd)


ggplot(sub_shannon, aes(x=mut_freq_pt_collapsed, y=Shannon_diversity)) + 
  geom_boxplot(notch=F, outlier.colour = "white", fill = "lightgrey") +
  geom_jitter(aes(fill = as.factor(high_sd_or_low_sd)), shape=21, position=position_jitter(0.2)) +
  scale_x_discrete(limits=c("2", "3", "4", "5", "6", "7" ,"8+")) +
  scale_fill_manual(values = c("navyblue", "darkred")) +
  labs(title = NULL) +
  ylab(label= "Shannon Diversity Index") +
  xlab(label = "Number of unique mutations per patient") +
  theme(legend.position="none") 

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/shannon_diversity_by_mutation_burden.pdf", dpi = 300, width = 5, height = 3, units = "in")



# make os time into years
sub_shannon$Time_to_OS <- (sub_shannon$Time_to_OS/365)

# create the survival data 
sub_shannon$OS <- with(sub_shannon, Surv(Time_to_OS, Censor == 1))

# create the survival objects used to plot kaplan-meyer curves
OS <- survfit(OS ~ high_sd_or_low_sd, data = sub_shannon, conf.type = "log-log")


# find the different p-values for the different comparisons
sub_shannon$Censor <- as.numeric(sub_shannon$Censor)  
res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ high_sd_or_low_sd,
                         data = sub_shannon)
print(res)

# plots the survival
surv_scatter_plot <- ggsurvplot(OS,
                                data = sub_shannon,
                                log = (OS),
                                log.rank.weights = c("survdiff"),
                                pval = T,
                                test.for.trend = F,
                                pval.method.size = 3,
                                pval.coord = c(0, 0),
                                conf.int = F,
                                censor = T,
                                surv.mean.line = "none",
                                risk.table = T,
                                risk.table.title = "",
                                risk.table.fontsize = 5,
                                risk.table.height = .25,
                                risk.table.y.text = T,
                                break.time.by = 1,
                                risk.table.pos = c("out"),
                                # palette = c("#a50f15", "#006d2c","#08519c", "darkgrey"),
                                title = "Survival by Shannon Diversity Index",
                                xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(20, "plain", "#252525"),
                                pval.size = 5,
                                font.x = c(20, "plain", "#252525"),
                                font.y =  c(20, "plain", "#252525"),
                                font.legend = c(15, "plain"),
                                font.tickslab = c(15, "plain", "#252525"),
                                # legend.labs = c("Clonal", "Subclonal RAS/RTK", "Subclonal CEBPA"),
                                legend.title = "Shannon\ndiversity quartile",
                                legend = "right",
                                ggtheme = theme(plot.title = element_text(hjust = 0.5, size = 15)))

print(surv_scatter_plot)

png(filename = "~/Desktop/survival_by_shannon_diversity.png", res = 300, width = 10, height = 7.5, units = "in")

surv_scatter_plot
print(surv_scatter_plot)
dev.off()

rm(list=setdiff(ls(), "final_data_matrix_2"))

