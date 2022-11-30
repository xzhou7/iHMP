no_function()
# set work directory

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

######work directory
masstools::setwd_project()
setwd("data_analysis/stool_microbiome_vs_metabolome/halla_example")

###next code should be run in the terminal in Rstudio
# cd /Users/xiaotaoshen/Box/Xiaotao Shen's Files/human_microbiome_project/data_analysis/stool_microbiome_vs_metabolome/halla_example
# halla -x data1.txt -y data2.txt -m spearman -o output_result

#######read result
result = read.table("output_result/all_associations.txt", header = TRUE)

#####use my method to calculate correlations
data1 = read.table("data1.txt", header = FALSE)
data2 = read.table("data2.txt", header = FALSE)

dim(data1)
dim(data2)

result2 = 
purrr::map(1:nrow(data1), function(idx1) {
  cat(idx1, " ")
  x = as.numeric(data1[idx1, -1])
  purrr::map(1:nrow(data1), function(idx2) {
    y = as.numeric(data2[idx2, -1])
    temp =
      cor.test(x, y, method = "spearman")
    data.frame(
      X_features = data1$V1[idx1],
      Y_features = data2$V1[idx2],
      association = unname(temp$estimate),
      p.values = unname(temp$p.value)
    )
  }) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

result2$q.values = p.adjust(result2$p.values, method = "BH")

dim(result)
dim(result2)
result$association == result2$association
plot(result$association, result2$association)
abline(0,1)
plot(result$p.values, result2$p.values)
abline(0,1)

plot(result$q.values, result2$q.values)
abline(0,1)


sum(result$q.values < 0.05)
sum(result2$q.values < 0.05)


