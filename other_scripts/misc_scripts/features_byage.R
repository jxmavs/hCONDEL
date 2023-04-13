library(tidyverse)
library(ggpubr)
library(magrittr)
library(rstatix)
library(reshape)
library(reshape2)
setwd("~/Desktop/MPRAdel project")

#plot 1: percent bar chart for feature broken down by age category
#read in condel table with deletion length, sequence MRCA, and feature
condel_table <- read.table("condel_agecollapse_feature.txt", sep="\t", header = TRUE, stringsAsFactors=TRUE)

condel_table$age <- factor(condel_table$age, levels = c("Euarchonta or younger", "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia", "Amniota", "Tetrapoda", "Vertebrata"))
condel_table$feature <- factor(condel_table$feature, levels = c("nongenic", "intronic", "UTR_3", "UTR_5", "TSS", "exon"))
table <- table(condel_table$age, condel_table$feature)
addmargins(table)
proportions <- round(prop.table(table,2),digits=4)*100
melt_proportions <- melt(proportions, id=age)

plot_percent_2 <- ggplot(melt_proportions, aes(x = Var2, y = value, fill = Var1)) + geom_col(position = "fill") + ylab("% Composition") + xlab("Feature Type") +
  scale_fill_manual(values=c("#777e9b", "#8ebddb", "#b3d5ea", "#d2e1ef", "#aeb6d9", "#c5b4d6", "#c3a0cc", "#ad84bb", "#967ea2"), name = "Age")
print(plot_percent_2)
ggsave("agepercent_byfeature.svg", plot=plot_percent_2)

#plot 2: replace age category with number for boxplot of deletion size (y) and age (x)
condel_table_2 <- read.table("condel_numeric_agecollapse_feature.txt", sep="\t", header = TRUE, stringsAsFactors=TRUE)
condel_table_2$age <- as.factor(condel_table_2$age)
plot_delsize_factor <- ggplot(condel_table_2, aes(x=age, y=del_len, fill = age)) + geom_jitter(alpha=0.2, width=0.3) + geom_boxplot(outlier.shape = NA, alpha = 0.8) + ylab("Deletion length, bp") + xlab("Time to MRCA (my)") +
  scale_fill_manual(values=c("#777e9b", "#8ebddb", "#b3d5ea", "#d2e1ef", "#aeb6d9", "#c5b4d6", "#c3a0cc", "#ad84bb", "#967ea2"), name = "Age", labels = c("Euarchonta or younger", "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia", "Amniota", "Tetrapoda", "Vertebrata"))
print(plot_delsize_factor)
ggsave("agevsdel_length.svg", plot=plot_delsize_factor)

pwc <- condel_table_2 %>% pairwise_t_test(del_len ~ age, pool.sd = FALSE)
pwc

#plot 3: genomic feature vs size, boxplot again (genomic feature x, size y)
condel_table_2$feature <- factor(condel_table_2$feature, levels = c("nongenic", "intronic", "UTR_3", "UTR_5", "TSS", "exon"))
pwc_2 <- condel_table_2 %>% pairwise_t_test(del_len ~ feature, pool.sd = FALSE)
pwc_2
pwc_2 <- pwc_2 %>% add_xy_position(x = "feature")
plot_delsize_age <- ggplot(condel_table_2, aes(x=feature, y=del_len)) + geom_jitter(alpha=0.2, width=0.3) + geom_boxplot(outlier.shape = NA, alpha = 0.8) + stat_pvalue_manual(pwc_2, hide.ns = TRUE) +
  ylab("Deletion length, bp") + xlab("Genomic Feature") 
print(plot_delsize_age)
ggsave("featurevsdel_length.svg", plot_delsize_age)
