
load("~/Desktop/MetaAML_results/final_data_matrix.RData")


survival_by_mut_burden_function <- function(pt_subset, mut_load_upper_cutoff, mut_load_lower_cutoff){
  
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix
    label1 <- "All"
  }
if(pt_subset == "De novo"){
  final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
  label1 <- "De novo"
}
if(pt_subset == "Secondary"){
  final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "transformed")
  label1 <- "Secondary"
}
if(pt_subset == "Relapse"){
  final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
  label1 <- "Relapse"
}
if(pt_subset == "Therapy related"){
  final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
  label1 <- "Therapy related"
}
if(pt_subset == "Other"){
  final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
  label1 <- "Other"
}
  
# add column for survival strate based on mutation burden
final_data_matrix_sub$strata <- NA


for(i in 1:nrow(final_data_matrix_sub)){
  if(final_data_matrix_sub$number_of_mutation[i] > mut_load_upper_cutoff){
    final_data_matrix_sub$strata[i] <- paste("Over", mut_load_upper_cutoff, "mutations")
  }
  if(final_data_matrix_sub$number_of_mutation[i] > mut_load_lower_cutoff & final_data_matrix_sub$number_of_mutation[i] < mut_load_upper_cutoff){
    final_data_matrix_sub$strata[i] <- paste("Between", mut_load_upper_cutoff, "and", mut_load_lower_cutoff)
  }
  if(final_data_matrix_sub$number_of_mutation[i] <= mut_load_lower_cutoff){
    final_data_matrix_sub$strata[i] <- paste("Under", mut_load_lower_cutoff, "mutations")
  }
}
  

final_data_matrix_sub$Time_to_OS <- (final_data_matrix_sub$Time_to_OS/365)

# create the survival data 
final_data_matrix_sub$OS <- with(final_data_matrix_sub, Surv(Time_to_OS, Censor == 1))

# create the survival objects used to plot kaplan-meyer curves
OS <- survfit(OS ~ strata, data = final_data_matrix_sub, conf.type = "log-log")


# find the different p-values for the different comparisons
final_data_matrix_sub$Censor <- as.numeric(final_data_matrix_sub$Censor)  
res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ strata,
                         data = final_data_matrix_sub)
print(res)


# plots the survival
surv_mut_burden <- ggsurvplot(OS,
                                data = final_data_matrix_sub,
                                log = (OS),
                                log.rank.weights = c("survdiff"),
                                pval = T,
                                test.for.trend = F,
                                pval.method.size = 3,
                                pval.coord = c(0, 0),
                                conf.int = F,
                                censor = T,
                                surv.median.line = "none",
                                risk.table = T,
                                risk.table.title = "",
                                risk.table.fontsize = 5,
                                risk.table.height = .25,
                                risk.table.y.text = F,
                                break.time.by = 1,
                                risk.table.pos = c("out"),
                                palette = c("#08519c", "#a50f15", "#006d2c"),
                                title = "MetaAML survival by mutation burden",
                                xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(20, "plain", "#252525"),
                                pval.size = 5,
                                font.x = c(20, "plain", "#252525"),
                                font.y =  c(20, "plain", "#252525"),
                                font.legend = c(15, "plain"),
                                font.tickslab = c(15, "plain", "#252525"),
                                # legend.labs = c(),
                                legend.title = "",
                                legend = "top",
                                ggtheme = theme(plot.title = element_text(hjust = 0.5)))

print(surv_mut_burden)
}

survival_by_mut_burden_function(pt_subset = "All", mut_load_lower_cutoff = 5, mut_load_upper_cutoff = 15)
