library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

load("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/Robject/Revision_MultiOmes_0509.RData")

num_subjectIDs <- n_distinct(clinic.df$SubjectID)
num_subjectIDs

# Count the number of times each SubjectID appears
id_counts <- table(clinic.df$SubjectID)

# Subset the data frame to only include rows where SubjectID appears at least 5 times
clinic.df_filtered <- clinic.df[clinic.df$SubjectID %in% names(id_counts)[id_counts >= 5], ]

table(clinic.df_filtered$SubjectID)
colnames(clinic.df_filtered)
library(dplyr)
library(tidyr)

# calculate correlation per SubjectID
cor_df <- clinic.df_filtered %>%
  group_by(SubjectID) %>%
  summarize(cor = cor(A1C, GLU, use = "complete.obs"),
            pval = cor.test(A1C, GLU)$p.value)

cor_df$p_adjusted <- p.adjust(cor_df$pval, method="BH")


ggplot(data = cor_df, aes(x = p_adjusted, y = cor, label = SubjectID)) +
  geom_point(size = 3) +
  labs(x = "p_adjusted", y = "cor") + ggtitle("correlation between A1C and FPG in each individual") + theme_cowplot()


# calculate mean A1C and mean GLU for each SubjectID
mean_df <- clinic.df_filtered %>% 
  group_by(SubjectID) %>% 
  summarise(mean_A1C = mean(A1C, na.rm = TRUE), mean_GLU = mean(GLU, na.rm = TRUE)) %>% 
  na.omit()

# calculate correlation between mean_A1C and mean_GLU
cor_df <- cor(mean_df$mean_A1C, mean_df$mean_GLU, method = "pearson", use = "pairwise.complete.obs")

# calculate p value
p_val <- cor.test(mean_df$mean_A1C, mean_df$mean_GLU, method = "pearson", use = "pairwise.complete.obs")$p.value

# create scatter plot with correlation coefficient and p value annotation

ggplot(mean_df, aes(x = mean_A1C, y = mean_GLU, label = SubjectID)) +
  geom_point() +
  geom_text(nudge_x = 0.1, nudge_y = 0.1) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", x = min(mean_df$mean_A1C), y = max(mean_df$mean_GLU), label = paste("cor =", round(cor_df, 2), ", p =", round(p_val, 2)), hjust = 0, vjust = 1) +
  labs(x = "Mean A1C", y = "Mean GLU", title = "Correlation between mean A1C and mean GLU by SubjectID") +
  theme_bw()

