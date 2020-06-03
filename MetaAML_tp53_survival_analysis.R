# MetaAML_tp53_survival_analysis
#
# Brooks Benard
# bbenard@stanford.edu
# 03/27/2020
#

### This script uses variant allele frequencies from MetaAML to understand clonal patterns at both genetic and survival levels

#### load required packages and data ####
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
if (!require('ggridges')) install.packages('ggridges'); library('ggridges')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('ggplotify')) install.packages('ggplotify'); library('ggplotify')
if (!require('survivalAnalysis')) install.packages('survivalAnalysis'); library('survivalAnalysis')
if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('vegan')) install.packages('vegan'); library('vegan')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('scales')) install.packages('scales'); library('scales')


# load the MetaAML results file 
load("~/Desktop/MetaAML_results/final_data_matrix.RData")
final_data_matrix_2 = subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
final_data_matrix_2 = na.omit(select(final_data_matrix_2, Sample, Gene, VAF, Censor, Time_to_OS, Subset))

# assign TP53 allelic status
for(i in 1:nrow(final_data_matrix_2)){
  if(final_data_matrix_2$Gene[i] == "TP53"){
    final_data_matrix_2$tp53_zygosity[i] = ifelse(final_data_matrix_2$VAF[i] >= .6, "loh", "mono-allelic")
  } else {
    final_data_matrix_2$tp53_zygosity[i] = "WT"
  }
}

# list of TP53 patients
tp53_pts = subset(final_data_matrix_2, final_data_matrix_2$Gene == "TP53")
tp53_pts = tp53_pts %>% group_by(Sample, Gene) %>% mutate(n_tp53 = n())
tp53_pts$tp53_zygosity = ifelse(tp53_pts$n_tp53 > 1, "bi-allelic", tp53_pts$tp53_zygosity)
tp53_pts = tp53_pts %>% select(Sample, tp53_zygosity, Time_to_OS, Censor) %>% unique()

# non-tp53 patients
wt_pts = subset(final_data_matrix_2, final_data_matrix_2$Gene != "TP53")
wt_pts = wt_pts %>% select(Sample, tp53_zygosity, Time_to_OS, Censor) %>% unique()

# combine
p53_all = rbind.fill(tp53_pts, wt_pts)

## Kaplan Meier analysis
# create the survival data object
p53_all$Time_to_OS <- (p53_all$Time_to_OS/365)

p53_all$OS <- with(p53_all, Surv(Time_to_OS, Censor == 1))

# create the survival objects used to plot kaplan-meyer curves
OS <- survfit(OS ~ tp53_zygosity, data = p53_all, conf.type = "log-log")

# find the different p-values for the different comparisons
p53_all$Censor = as.numeric(p53_all$Censor)
res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ tp53_zygosity,
                         data = p53_all)
print(res)

# cohorts <- c("", "2", "3", "4")

# plots the survival
surv_plot <- ggsurvplot(OS,
                        data = p53_all,
                        log = (OS),
                        log.rank.weights = c("survdiff"),
                        pval = F,
                        test.for.trend = F,
                        pval.method.size = 3,
                        pval.coord = c(0, 0),
                        conf.int = F,
                        censor = T,
                        surv.median.line =  "hv",
                        risk.table = T,
                        risk.table.title = "",
                        risk.table.fontsize = 5,
                        risk.table.height = .3,
                        risk.table.y.text = T,
                        break.time.by = 1,
                        risk.table.pos = c("out"),
                        palette = c("#0570b0", "#238b45", "#cb181d","#bdbdbd"),
                        xlab = "Years",
                        ylim = c(0, 1.0),
                        xlim = c(0,5),
                        ylab =  "Survival Probability",
                        font.main = c(15, "plain", "#252525"),
                        pval.size = 5,
                        font.x = c(15, "plain", "#252525"),
                        font.y =  c(15, "plain", "#252525"),
                        font.legend = c(15, "plain"),
                        font.tickslab = c(15, "plain", "#252525"),
                        legend.labs = c("bi-allelic", "loh", "mono-allelic", "WT"),
                        legend.title = "TP53 status",
                        legend = "right",
                        ggtheme = theme(plot.title = element_text(hjust = 0.5)))

print(surv_plot)

 png(filename = "~/Desktop/MetaAML_results/Data/Figures/meta_aml_survival_by_TP53_detailed.png", res = 300, width = 8, height = 7.5, units = "in")

surv_plot
print(surv_plot)
dev.off()
