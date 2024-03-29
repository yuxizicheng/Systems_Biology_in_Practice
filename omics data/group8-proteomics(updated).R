proteomics_data <- read.csv("/Users/liucheng/Desktop/omics data/Proteomics (1).csv", header = TRUE)
# remove rows with "ratio_count=0"
proteomics_data = subset(proteomics_data, ratio_count != 0)
proteomics_data
# remove missing values
proteomics_data <- na.omit(proteomics_data)

# calculate the t-value for each protein
proteomics_data$t_value = proteomics_data$avg_ratio/(proteomics_data$sd_ratio/sqrt(proteomics_data$ratio_count))
# Calculate the p-value for each protein using a two-sided t-test
proteomics_data$p_value = 2* pt(abs(proteomics_data$t_value), df = proteomics_data$ratio_count-1, lower.tail = FALSE)
# adjust the p-values for multiple comparisons
proteomics_data$p_adjusted = p.adjust(proteomics_data$p_value, method = "fdr")
# view the results
print(proteomics_data[,c("protein", "p_value", "p_adjusted")])
proteomics_data

library(ggplot2)
#install.packages("ggrepel")
library(ggrepel) 
proteomics_data$label=ifelse(abs(proteomics_data$avg_ratio)>1.5 & proteomics_data$p_adjusted<0.1 , as.character(proteomics_data$protein), '')
# create a data frame with fold change and p_adjusted columns
volcano_data <- data.frame(log2_fc =proteomics_data$avg_ratio, neg_log10pvalue = -log10(proteomics_data$p_adjusted),protein=proteomics_data$protein,label=proteomics_data$label)
# create a volcano plot
volcano_plot = ggplot(data = volcano_data, aes(x= log2_fc, y=neg_log10pvalue)) +
  geom_point(alpha = 0.5, size=2) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black",size=0.7) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black" ,size=0.7) +
  labs(title = "The Volcano Plot of Proteomics", x = "log2 Fold Change", y ="-log10(Adjusted P-value)") +
  theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16),
        legend.title = element_blank()) +
  
  scale_x_continuous(limits = c(-5,5)) 

# add significant proteins as a new layer
  sig_proteins = proteomics_data[abs(proteomics_data$avg_ratio)>1.5 &proteomics_data$p_adjusted<0.1,]
  sig_proteins_data = data.frame(log2_fc = sig_proteins$avg_ratio, neg_log10pvalue = -log10(sig_proteins$p_adjusted),protein=sig_proteins$protein)
  sig_proteins_layer = geom_point(data = sig_proteins_data, aes(x=log2_fc, y=neg_log10pvalue),color= "red",size = 3, alpha = 0.8)
  
  
volcano_plot_with_sig_proteins =volcano_plot +sig_proteins_layer +
  scale_color_manual(values = c("Differentially Expressed Proteins" = "red")) +
  guides(color = guide_legend(title = "Legend Title"))+
  geom_text_repel(aes(x = log2_fc,                   # geom_text_repel 标记函数
                      y = neg_log10pvalue,          
                      label=volcano_data$label),                       
                  max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
                  size=3,                                  # 字体大小
                  box.padding=unit(0.5,'lines'),           # 标记的边距
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                   # 标记线条的颜色
                  show.legend=TRUE)
volcano_plot_with_sig_proteins

sig_proteins_info= as.data.frame(sig_proteins$protein,sig_proteins$description)
sig_proteins_info



























