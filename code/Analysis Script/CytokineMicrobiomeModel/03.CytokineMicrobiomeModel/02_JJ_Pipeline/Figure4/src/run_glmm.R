# 01_model.R
# 
# 2018-09-30
# Jethro (modified from DS code)
# 

suppressMessages(library('tidyverse'))
suppressMessages(library('brms'))
suppressMessages(library('broom'))

## Define input variables. 
# 1) model_table
# 2) tax_table
# 3) name of desired taxon
# 4) name of desired cytokine
# 5) output file name
# 6) file in which to store full brm model

args <- commandArgs(TRUE)
model_df = args[1]
tax_df = args[2]
taxon = args[3]
cytokine = args[4]
outfile = args[5]
modfile = args[6]

# print(model_df)
# print(tax_df)
# print(taxon)
# print(cytokine)
# print(outfile)

### Load the existing R data
# model.df
load(model_df)
# df.tax
load(tax_df)

### Subset taxon table so that it contains only taxon of interest
tax.df = df.tax[, names(df.tax) %in% c(taxon, 'SampleID')]
model.df = merge(model.df, tax.df, by='SampleID')

# Filter out cluster 2
x <- model.df %>%
  filter(cluster != 2) %>%
  mutate(cluster = fct_drop(cluster))

print(head(x))

### Run model according to DS code
c <- cytokine
g <- taxon

f <- as.formula(paste(g, "~ days + (1|SubjectID) + cluster *", c, sep = " "))
    
brm <- brms::brm(f, data = x, sparse = TRUE, family = negbinomial())

saveRDS(brm, file = modfile)


# Reformat output
brm.df <- broom::tidy(brm)
    
# Create column to track the genus being modelled
brm.df$genus <- g

# Create list of outputs for each genus
model.output <- 
      brm.df %>%
      mutate(significant = !lower <= 0 & upper >= 0) %>%
      filter(grepl(c, brm.df$term))

saveRDS(model.output, file = outfile)