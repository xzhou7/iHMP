####
#combineing the plot 
library(ggplot2)
library(dplyr)

setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/human_microbiome_project")
getwd()

load("./data_analysis/proteome/data_preparation/cor_data")
cor_data_pro <- cor_data

plot <- cor_data_pro %>% 
  ggplot(aes(cor)) +
  geom_histogram(color = "black", binwidth = 0.05) +
  theme_bw() +
  labs(x = "Spearman correlation",
       y = "Count") + ggtitle("Proteome")
plot
ggsave(filename = "./Figures/FigureS1/protein_correlation distributation.pdf", plot, height = 5, width = 5, dpi = 300)

load("./data_analysis/metabolome/data_preparation/cor_data")
cor_data_met <- cor_data

plot2 <- cor_data_met %>% 
  ggplot(aes(cor)) +
  geom_histogram(color = "black", binwidth = 0.05) +
  theme_bw() +
  labs(x = "Spearman correlation",
       y = "Count") + ggtitle("Metabolome")
plot2
ggsave(filename = "./Figures/FigureS1/metabolite_correlation distributation.pdf", plot2, height = 5, width = 5, dpi = 300)

load("./data_analysis/lipidome/data_preparation/cor_data")
cor_data_lip <- cor_data

plot3 <- cor_data_lip %>% 
  ggplot(aes(cor)) +
  geom_histogram(color = "black", binwidth = 0.05) +
  theme_bw() +
  labs(x = "Spearman correlation",
       y = "Count") + ggtitle("Lipidome")
plot3
ggsave(filename = "./Figures/FigureS1/lipid_correlation distributation.pdf", plot3, height = 5, width = 5, dpi = 300)






