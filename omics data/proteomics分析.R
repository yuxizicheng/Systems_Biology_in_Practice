if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)
install.packages("caret")
install.packages("vctrs")
library(caret)
# load data from file
proteomics_data <- read.csv("/Users/liucheng/Desktop/omics data/Proteomics (1).csv", header = TRUE)
# remove missing values
proteomics_data <- na.omit(proteomics_data)
# normalize data
proteomics_data_norm = (proteomics_data[,3:5]-proteomics_data[,6])/proteomics_data[,7]
# summary statistics
summary(proteomics_data_norm)
# histogram of log_ratio_1
hist(proteomics_data_norm[,3])
proteomics_data_norm









boxplot(proteomics_data$log_ratio_1) 
hist(proteomics_data$log_ratio_3)
boxplot(proteomics_data$avg_ratio, horizontal= TRUE)
plot(proteomics_data$log_ratio_1, proteomics_data$log_ratio_2)
cor(proteomics_data$log_ratio_1, proteomics_data$log_ratio_2) # 0.6518216
cor(proteomics_data$log_ratio_1, proteomics_data$log_ratio_3) #0.664428
cor(proteomics_data$log_ratio_2, proteomics_data$log_ratio_3) #0.7825619
cor_matrix<-cor(proteomics_data[,c('log_ratio_1','log_ratio_2','log_ratio_3')])
highly_correlated = findCorrelation(cor_matrix, cutoff = 0.8)

# Pre-process the data
proteomics_data <- na.omit(proteomics_data)
proteomics_data$log_ratio_1 <- log2(proteomics_data$log_ratio_1) # transform the data if necessary
proteomics_data$log_ratio_2 <- log2(proteomics_data$log_ratio_2)
proteomics_data$log_ratio_3 <- log2(proteomics_data$log_ratio_3)

# Perform statistical analysis
proteomics_results <- t.test(log_ratio_1 ~ log_ratio_3, data = proteomics_data) # perform t-test
proteomics_results <- data.frame(proteins = rownames(proteomics_results$estimate), # create a data frame of results
                                 p_value = proteomics_results$p.value,
                                 avg_ratio = proteomics_data$avg_ratio,
                                 log_ratio_1 = proteomics_data$log_ratio_1,
                                 log_ratio_2 = proteomics_data$log_ratio_2,
                                 log_ratio_3 = proteomics_data$log_ratio_3,
                                 ratio_count = proteomics_data$ratio_count)

# Correct for multiple testing
proteomics_results$p_adjusted <- p.adjust(proteomics_results$p_value, method = "bonferroni") # perform Bonferroni correction
proteomics_results$p_adjusted_fdr <- p.adjust(proteomics_results$p_value, method = "fdr") # perform FDR correction


# Visualize the results
library(ggplot2)
ggplot(proteomics_results, aes(x = log_ratio_2, y = log_ratio_3, color = p_adjusted_fdr)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  xlab("Log Ratio 2") +
  ylab("Log Ratio 3") +
  ggtitle("Volcano Plot of Proteomics Data")

























# convert log ratios to matrix and set column names
log_ratios = data.matrix(proteomics_data_norm[,c("log_ratio_1","log_ratio_2","log_ratio_3")])
colnames(log_ratios) = c("A","B","C")

# create design matrix including all variables
#design = model.matrix(~0 + cbind(log_ratios), data = proteomics_data_norm, weights= proteomics_data$ratio_count)
design <- model.matrix(~ 0 + colnames(log_ratios))


# fit linear model
fit = lmFit(log_ratios, design)

# check for missing values in the response variable
sum(is.na(fit$Y)) # count the number of missing values

# calculate the residual standard deviation
summary(fit)$sigma


# estimate variance using empirical Bayes
fit = eBayes(fit)

# identify differentially expressed proteins based on all variables
results = topTable(fit, coef=1, number=50, sort.by='p', p.value= 0.05)

#view the top differentially expressed proteins
head(results)

# proteomics_data <- read.csv("/Users/liucheng/Desktop/omics data/Proteomics (1).csv", header = TRUE)
# library(ggplot2)
# ggplot(proteomics_data, aes(x = avg_ratio, y = sd_ratio)) +
#   geom_point() +
#   labs(x = "Log2 fold change", y = "Standard deviation", title = "Volcano plot of differential proteins")
# threshold <- with(proteomics_data, abs(avg_ratio) > 1.5 & ratio_count == 3 & pvalue < 0.05)
# ggplot(proteomics_data, aes(x = avg_ratio, y = sd_ratio)) +
#   geom_point(aes(color = threshold)) +
#   labs(x = "Log2 fold change", y = "Standard deviation", title = "Volcano plot of differential proteins") +
#   scale_color_manual(values = c("gray", "red"), guide = FALSE) +
#   geom_hline(yintercept = 1.5, linetype = "dashed") +
#   geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed")
# 
# # Define the threshold for differential expression
# threshold <- with(proteomics_data, abs(avg_ratio) > 1.5 & ratio_count == 3 & pvalue < 0.05)
# 
# # Create the volcano plot
# ggplot(proteomics_data, aes(x = avg_ratio, y = sd_ratio)) +
#   geom_point(aes(color = threshold)) +
#   labs(x = "Log2 fold change", y = "Standard deviation", title = "Volcano plot of differential proteins") +
#   scale_color_manual(values = c("gray", "red"), guide = FALSE) +
#   geom_hline(yintercept = 1.5, linetype = "dashed") +
#   geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +
# 
#   # Add text labels for differentially expressed proteins
#   geom_text(aes(label = ifelse(threshold, protein, "")), size = 3, vjust = -0.5, hjust = 0.5)



