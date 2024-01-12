library(survival)
library(survminer)

#To set working directory on windows machine, replace backslashes with double backslashes for R to accept our file path
wd <- r"(C:\Users\mattm\OneDrive\Desktop\GitHub_projects\Survival_analyses_with_R)"

#set our working directory (where our files are stored and outputted)
setwd(wd)

#read in survival dataframe
#survival dataframe has three 
survival_df <- read.delim("sgArid1a_B16F10_isotype_vs_ICB.txt", header=T, sep = "\t" )

#check if group column is a factor
is.factor(survival_df$group)
#FALSE

#Make group column values a factor so we can put the variables in the order we prefer on our plot
survival_df$group <- as.factor(survival_df$group)
levels(survival_df$group)
#[1] "ICB"     "Isotype"

#Reorder factor level
survival_df$group <- factor(survival_df$group, levels = c("Isotype", "ICB"))
levels(survival_df$group)


#Fit survival model using survfit() function from survival package using time (v14) and censored (v15) and specify that values of 1 are events (death).
#We also specify that we want to pull our data from the sample dataset and parse the data by values in the ARID1A column.
fit <- survfit(Surv(time, as.numeric(substr(censor, 1, 1))) ~ group,
               data = survival_df
)
plot(fit)

#perform log-rank survival test (pariwise survival test)
log_rank_test <- survdiff(Surv(time, as.numeric(substr(censor, 1, 1))) ~ group, data = survival_df)
#output the p value from the log rank test
p_value <- round(1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1), 5)


#To make our kaplan meier curve plot more informative, let's add a color code, label the axes, a table below the plot to show the number of patients across time, a p-value from log-rank test, and confidence intervals. 
#This should make for a much more informative kaplan meier plot! 
survplot <- ggsurvplot(fit, data = survival_df, palette = c("red", "blue"), title = "sgArid1a B16F10 survival +/- Immunotherapy",
                       legend.title="Tumor Genotype", legend.labs=c("sgArid1a + ICB", "sgArid1a + Isotype Control"), 
                       pval = T, risk.table = TRUE, xlab = "Days Post Injection", conf.int = F)
survplot