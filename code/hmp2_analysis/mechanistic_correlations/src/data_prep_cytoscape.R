
library(dplyr)
library(reshape2)

# You can effectively remove scientific notation in printing with this code:
options(scipen=999)

# TODO: Change paths
D1_path <-  '~/Dropbox (Huttenhower Lab)/Faroese/Analysis/halla_output/halla_time_28_species/Distance_Matrix1.tsv' #~/Documents/Hutlab/halla/halla_output/
D2_path <-  '~/Dropbox (Huttenhower Lab)/Faroese/Analysis/halla_output/halla_time_28_species/Distance_Matrix2.tsv'
sim_path <- '~/Dropbox (Huttenhower Lab)/Faroese/Analysis/halla_output/halla_time_28_species/similarity_table.txt'

# can set these filters in cytoscape
within_corr_threshold <- .75
between_corr_threshold <- .29

# dm1
D1 <- read.table(D1_path, header = T,
                        row.names =1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
D1 <- 1- D1 # spearman
D1[lower.tri(D1, diag = TRUE)] <-NA
D1_data <- melt(D1 %>% mutate(source = rownames(.)), id.vars = 'source')
D1_data <- setNames(D1_data, c('source', 'destination', 'correlation'))
D1_data <- D1_data[abs(D1_data$correlation) > within_corr_threshold,]

D2 <- read.table(D2_path, header = T,
                 row.names =1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
D2 <- 1- D2
D2[lower.tri(D2, diag = TRUE)] <-NA
D2_data <- melt(D2 %>% mutate(source = rownames(.)), id.vars = 'source')
D2_data <- setNames(D2_data, c('source', 'destination', 'correlation'))

D2_data <- D2_data[abs(D2_data$correlation) > within_corr_threshold,]

sim_table <- read.table(sim_path, header = T,
                       row.names =1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)

sim_table_data <- melt(sim_table %>% mutate(source = rownames(.)), id.vars = 'source')
sim_table_data <- setNames(sim_table_data, c('source', 'destination', 'correlation'))
sim_table_data <- rbind(sim_table_data, D1_data)
sim_table_data <- rbind(sim_table_data, D2_data)

sim_table_data$correlation <- as.numeric(sim_table_data$correlation)
sim_table_data <- sim_table_data[abs(sim_table_data$correlation) > between_corr_threshold,]

sim_table_data <- sim_table_data[sim_table_data$source != sim_table_data$destination,]
#sim_table_data <- sim_table_data[complete.cases(sim_table_data),]
sim_table_data <- sim_table_data[!is.na(sim_table_data$correlation),]
sim_table_data$color <- NA
sim_table_data$color[sim_table_data$correlation>0] <- 'Red'
sim_table_data$color[sim_table_data$correlation<=0] <- 'Blue'
sim_table_data$type <- 'cooccurance'
sim_table_data$correlation <- abs(sim_table_data$correlation) * 10



#sim_table_data <- sim_table_data[c('source', 'person', 'time', 'Y')]
write.table(sim_table_data, '~/Documents/Faroese_project/halla_output/cytoscape_data.txt', sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
