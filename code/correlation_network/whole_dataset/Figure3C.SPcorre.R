correlation.df <- read.csv("/Users/xzhou7/Library/CloudStorage/Box-Box/human_microbiome_project/Figures/Figure3/detailed/figure3c_data/beta.estimation.fourbodysite.csv", header = T)


cor.test(correlation.df$stool, correlation.df$skin, method = "spearman")
cor.test(correlation.df$stool, correlation.df$oral, method = "spearman")
cor.test(correlation.df$stool, correlation.df$nasal, method = "spearman")

cor.test(correlation.df$skin, correlation.df$oral, method = "spearman")
cor.test(correlation.df$skin, correlation.df$nasal, method = "spearman")

cor.test(correlation.df$oral, correlation.df$nasal, method = "spearman")
