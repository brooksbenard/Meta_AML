# MetaAML_cohort_summary
# Brooks Benard
# bbenard@stanford.edu
# 04.22.2020
# this script visualized indicidual cohort level aspects such as age, survival, mutations, vafs, etc. for the Meta-AML dataset

# packages
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('ggridges')) install.packages('ggridges'); library('ggridges')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')

# data
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

# establish uniform colors for each cohort
cohort_colors = c("Tyner" = "#0073C2FF", 
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

# Age ####
sub=na.omit(as.data.frame(distinct(final_data_matrix, Sample, Age, Cohort)))
sub$Age=as.numeric(sub$Age)

ggplot(sub, aes(x = Age, y = as.factor(Cohort), fill = Cohort)) +
  geom_density_ridges(
    jittered_points = F, position = "raincloud",
    alpha = 1, scale = 1,
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_y_discrete(expand = c(0,0)) +
  ylab(label = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values = cohort_colors)

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/age_by_cohort.png", dpi = 300, width = 4, height = 5, units = "in")

# Gender ####
sub=na.omit(as.data.frame(distinct(final_data_matrix, Sample, Sex, Cohort)))

sub = sub %>% group_by(Cohort, Sex) %>% tally()

ggplot(sub, aes(fill=Sex, y=n, x=Cohort)) + 
  geom_bar(position = "fill", stat="identity",) +
  scale_fill_manual(values = c("Male" = "#6a51a3", "Female" = "#43a2ca")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  ylab(label = NULL) +
  xlab(label = NULL) 
  
ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/sex_by_cohort.png", dpi = 300, width = 5, height = 3, units = "in")

# Risk ####
sub=na.omit(as.data.frame(distinct(final_data_matrix, Sample, Risk, Cohort)))

sub = sub %>% group_by(Cohort, Risk) %>% tally()

ggplot(sub, aes(fill=Risk, y=n, x=Cohort)) + 
  geom_bar(position = "fill", stat="identity",) +
  scale_fill_manual(values = c("Adverse" = "#E64B35FF",
                               "Intermediate" = "#8491B4FF",
                               "Favorable" = "#00A087FF", 
                               "Unknown" = "#767676FF")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   # size = 12, 
                                   hjust = 1)) +
  ylab(label = NULL) +
  xlab(label = NULL) 

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/risk_by_cohort.png", dpi = 300, width = 5, height = 3, units = "in")


 # Survival ####
sub = na.omit(as.data.frame(distinct(final_data_matrix, Sample, Cohort, Time_to_OS, Censor)))
sub$Cohort = as.factor(sub$Cohort)

# all
sub$Time_to_OS <- (sub$Time_to_OS/365)
sub$OS <- with(sub, Surv(Time_to_OS, Censor == 1))
OS <- survfit(OS ~ Cohort, data = sub, conf.type = "log-log")

surv_plot <- ggsurvplot(OS,
                        data = sub,
                        log = (OS),
                        log.rank.weights = c("survdiff"),
                        pval = F,
                        test.for.trend = F,
                        pval.method.size = 3,
                        pval.coord = c(0, 0),
                        conf.int = F,
                        censor = T,
                        surv.median.line = "none",
                        risk.table = F,
                        risk.table.title = "",
                        risk.table.fontsize = 4,
                        risk.table.height = .3,
                        risk.table.y.text = T,
                        break.time.by = 5,
                        risk.table.pos = c("out"),
                        palette = cohort_colors,
                        xlab = "Years",
                        ylim = c(0, 1.0),
                        ylab =  "Survival Probability",
                        font.main = c(15, "plain", "#252525"),
                        pval.size = 4,
                        font.x = c(12, "plain", "#252525"),
                        font.y =  c(12, "plain", "#252525"),
                        font.legend = c(12, "plain"),
                        font.tickslab = c(12, "plain", "#252525"),
                        legend.labs = c("Garg", "Hirsch", "Huet", "Lindsley", "Majeti", "Papaemmanuil", "TCGA", "Tyner", "Welch"),
                        legend.title = "Cohort",
                        legend = "right",
                        # title = gene,
                        ggtheme = theme(plot.title = element_text(hjust = 0.5)))

print(surv_plot)
png(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/total_survival.png", res = 300, width = 7.5, height = 5, units = "in")

surv_plot
print(surv_plot)
dev.off()

# de novo
sub = na.omit(as.data.frame(distinct(final_data_matrix, Sample, Cohort, Time_to_OS, Censor, Subset)))

sub = subset(sub, Subset == "de_novo")

sub$Time_to_OS <- (sub$Time_to_OS/365)
sub$OS <- with(sub, Surv(Time_to_OS, Censor == 1))
OS <- survfit(OS ~ Cohort, data = sub, conf.type = "log-log")

surv_plot <- ggsurvplot(OS,
                        data = sub,
                        log = (OS),
                        log.rank.weights = c("survdiff"),
                        pval = F,
                        test.for.trend = F,
                        pval.method.size = 3,
                        pval.coord = c(0, 0),
                        conf.int = F,
                        censor = T,
                        surv.median.line = "none",
                        risk.table = F,
                        risk.table.title = "",
                        risk.table.fontsize = 4,
                        risk.table.height = .3,
                        risk.table.y.text = T,
                        break.time.by = 5,
                        risk.table.pos = c("out"),
                        palette = cohort_colors,
                        xlab = "Years",
                        ylim = c(0, 1.0),
                        ylab =  "Survival Probability",
                        font.main = c(15, "plain", "#252525"),
                        pval.size = 4,
                        font.x = c(12, "plain", "#252525"),
                        font.y =  c(12, "plain", "#252525"),
                        font.legend = c(12, "plain"),
                        font.tickslab = c(12, "plain", "#252525"),
                        legend.labs = c("Garg", "Hirsch", "Huet", "Majeti", "Papaemmanuil", "TCGA", "Tyner", "Welch"),
                        legend.title = "Cohort",
                        legend = "right",
                        # title = gene,
                        ggtheme = theme(plot.title = element_text(hjust = 0.5)))

print(surv_plot)
png(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/de_novo_survival.png", res = 300, width = 7.5, height = 5, units = "in")

surv_plot
print(surv_plot)
dev.off()


# Number of mutations per patient ####
# plot total, variant type, and risk



# Most frequent mutations ####
# plot total, variant type, and risk
sub = na.omit(as.data.frame(distinct(final_data_matrix, Sample, Gene, VAF, Cohort, Subset)))


# VAF per cohort ####
sub = na.omit(as.data.frame(distinct(final_data_matrix, Sample, Gene, VAF, Cohort, Subset)))

ggplot(sub, aes(x = VAF, y = as.factor(Cohort), fill = Cohort)) +
  geom_density_ridges(
    jittered_points = F, 
    position = "raincloud",
    alpha = 1, scale = 1,
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_y_discrete(expand = c(0,0)) +
  ylab(label = NULL) +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values = cohort_colors)

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/vaf_by_cohort.png", dpi = 300, width = 4, height = 5, units = "in")


# VAF per gene ####
sub <- subset(final_data_matrix, final_data_matrix$mut_freq_gene >= 100)
sub = na.omit(as.data.frame(distinct(sub, Sample, Gene, VAF, Cohort, Subset)))

ggplot(sub, aes(x = VAF, y = reorder(Gene, desc(Gene)))) +
  geom_density_ridges(
    jittered_points = F, 
    position = "raincloud",
    alpha = 1, scale = 1,
    quantile_lines = TRUE, quantiles = 2
  ) +
  scale_y_discrete(expand = c(0,0)) +
  xlim(0,1) +
  ylab(label = NULL) +
  theme(legend.position = "none") 

ggsave(filename = "~/Desktop/MetaAML_results/Data/Figures/cohort_summary/vaf_by_gene.png", dpi = 300, width = 4, height = 5, units = "in")
