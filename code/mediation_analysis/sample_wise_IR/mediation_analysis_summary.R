####reference 
####https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171

# no_function()
# set work directory

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

load("data_analysis/stool_microbiome/data_preparation/variable_info")
stool_microbiome_variable_info <- variable_info

load("data_analysis/skin_microbiome/data_preparation/variable_info")
skin_microbiome_variable_info <- variable_info

load("data_analysis/oral_microbiome/data_preparation/variable_info")
oral_microbiome_variable_info <- variable_info

load("data_analysis/nasal_microbiome/data_preparation/variable_info")
nasal_microbiome_variable_info <- variable_info

load("data_analysis/proteome/data_preparation/variable_info")
proteome_variable_info <- 
  variable_info %>% 
  dplyr::select(variable_id) %>% 
  dplyr::mutate(true_name = variable_id)

load("data_analysis/metabolome/data_preparation/variable_info")
metabolome_variable_info <- 
  variable_info %>% 
  dplyr::select(variable_id, Metabolite) %>% 
  dplyr::rename(true_name = Metabolite)

load("data_analysis/cytokine/data_preparation/variable_info")
cytokine_variable_info <- 
  variable_info %>% 
  dplyr::select(variable_id) %>% 
  dplyr::mutate(true_name = variable_id)

load("data_analysis/lipidome/data_preparation/new_variable_info")
lipidome_variable_info <- 
  new_variable_info %>% 
  dplyr::select(variable_id, annotation) %>% 
  dplyr::rename(true_name = annotation)

####stool microbiome
####cytokine
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_cytokine/mediation_result")
mediation_result <-
  mediation_result %>% 
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_cytokine/mediation_result_inverse")

mediation_result <- 
mediation_result %>% 
  dplyr::left_join(stool_microbiome_variable_info, by = c("treat" = "variable_id")) %>% 
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <- 
mediation_result %>% 
  dplyr::left_join(mediation_result_inverse, 
                   by = c("phenotype", "mediator", "treat"))

stool_microbiome_cytokine_mediation_result <- 
  mediation_result

stool_microbiome_cytokine_mediation_result <- 
stool_microbiome_cytokine_mediation_result %>% 
  dplyr::mutate(treat_class = "Stool microbiome",
                mediator_class = "Cytokine",
                phenotype_class = "Phenotype") %>% 
  dplyr::left_join(cytokine_variable_info, 
                   by = c("mediator" = "variable_id")) %>% 
  dplyr::rename(mediator_true_name = true_name) %>% 
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####proteome
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_proteome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_proteome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(stool_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

stool_microbiome_proteome_mediation_result <-
  mediation_result

stool_microbiome_proteome_mediation_result <-
  stool_microbiome_proteome_mediation_result %>%
  dplyr::mutate(treat_class = "Stool microbiome",
                mediator_class = "Proteome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(proteome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####lipidome
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_lipidome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05 & prop_mediate < 1)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_lipidome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(stool_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

stool_microbiome_lipidome_mediation_result <-
  mediation_result

stool_microbiome_lipidome_mediation_result <-
  stool_microbiome_lipidome_mediation_result %>%
  dplyr::mutate(treat_class = "Stool microbiome",
                mediator_class = "Lipidome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(lipidome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####metabolome
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_metabolome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/stool_microbiome_vs_metabolome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(stool_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

stool_microbiome_metabolome_mediation_result <-
  mediation_result

stool_microbiome_metabolome_mediation_result <-
  stool_microbiome_metabolome_mediation_result %>%
  dplyr::mutate(treat_class = "Stool microbiome",
                mediator_class = "Metabolome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(metabolome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)


####skin microbiome
####cytokine
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_cytokine/mediation_result")
mediation_result <-
  mediation_result %>% 
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_cytokine/mediation_result_inverse")

mediation_result <- 
  mediation_result %>% 
  dplyr::left_join(skin_microbiome_variable_info, by = c("treat" = "variable_id")) %>% 
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <- 
  mediation_result %>% 
  dplyr::left_join(mediation_result_inverse, 
                   by = c("phenotype", "mediator", "treat"))

skin_microbiome_cytokine_mediation_result <- 
  mediation_result

skin_microbiome_cytokine_mediation_result <- 
  skin_microbiome_cytokine_mediation_result %>% 
  dplyr::mutate(treat_class = "Skin microbiome",
                mediator_class = "Cytokine",
                phenotype_class = "Phenotype") %>% 
  dplyr::left_join(cytokine_variable_info, by = c("mediator" = "variable_id")) %>% 
  dplyr::rename(mediator_true_name = true_name) %>% 
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####proteome
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_proteome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_proteome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(skin_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

skin_microbiome_proteome_mediation_result <-
  mediation_result

skin_microbiome_proteome_mediation_result <-
  skin_microbiome_proteome_mediation_result %>%
  dplyr::mutate(treat_class = "Skin microbiome",
                mediator_class = "Proteome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(proteome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)


####lipidome
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_lipidome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_lipidome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(skin_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

skin_microbiome_lipidome_mediation_result <-
  mediation_result


skin_microbiome_lipidome_mediation_result <-
  skin_microbiome_lipidome_mediation_result %>%
  dplyr::mutate(treat_class = "Skin microbiome",
                mediator_class = "Lipidome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(lipidome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)


####metabolome
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_metabolome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/skin_microbiome_vs_metabolome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(skin_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

skin_microbiome_metabolome_mediation_result <-
  mediation_result

skin_microbiome_metabolome_mediation_result <-
  skin_microbiome_metabolome_mediation_result %>%
  dplyr::mutate(treat_class = "Skin microbiome",
                mediator_class = "Metabolome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(metabolome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####oral microbiome
####cytokine
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_cytokine/mediation_result")
mediation_result <-
  mediation_result %>% 
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_cytokine/mediation_result_inverse")

mediation_result <- 
  mediation_result %>% 
  dplyr::left_join(oral_microbiome_variable_info, by = c("treat" = "variable_id")) %>% 
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <- 
  mediation_result %>% 
  dplyr::left_join(mediation_result_inverse, 
                   by = c("phenotype", "mediator", "treat"))

oral_microbiome_cytokine_mediation_result <- 
  mediation_result

oral_microbiome_cytokine_mediation_result <- 
  oral_microbiome_cytokine_mediation_result %>% 
  dplyr::mutate(treat_class = "Oral microbiome",
                mediator_class = "Cytokine",
                phenotype_class = "Phenotype") %>% 
  dplyr::left_join(cytokine_variable_info, by = c("mediator" = "variable_id")) %>% 
  dplyr::rename(mediator_true_name = true_name) %>% 
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####proteome
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_proteome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_proteome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(oral_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

oral_microbiome_proteome_mediation_result <-
  mediation_result

oral_microbiome_proteome_mediation_result <-
  oral_microbiome_proteome_mediation_result %>%
  dplyr::mutate(treat_class = "Oral microbiome",
                mediator_class = "Proteome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(proteome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)


####lipidome
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_lipidome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_lipidome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(oral_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

oral_microbiome_lipidome_mediation_result <-
  mediation_result

oral_microbiome_lipidome_mediation_result <-
  oral_microbiome_lipidome_mediation_result %>%
  dplyr::mutate(treat_class = "Oral microbiome",
                mediator_class = "Lipidome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(lipidome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####metabolome
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_metabolome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/oral_microbiome_vs_metabolome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(oral_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

oral_microbiome_metabolome_mediation_result <-
  mediation_result

oral_microbiome_metabolome_mediation_result <-
  oral_microbiome_metabolome_mediation_result %>%
  dplyr::mutate(treat_class = "Oral microbiome",
                mediator_class = "Metabolome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(metabolome_variable_info,
                   by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####nasal microbiome
####cytokine
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_cytokine/mediation_result")
mediation_result <-
  mediation_result %>% 
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_cytokine/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(nasal_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

nasal_microbiome_cytokine_mediation_result <-
  mediation_result

nasal_microbiome_cytokine_mediation_result <-
  nasal_microbiome_cytokine_mediation_result %>%
  dplyr::mutate(treat_class = "Nasal microbiome",
                mediator_class = "Cytokine",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(cytokine_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
dplyr::mutate(treat_true_name = Genus,
              phenotype_true_name = phenotype)

####proteome
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_proteome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_proteome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(nasal_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

nasal_microbiome_proteome_mediation_result <-
  mediation_result

nasal_microbiome_proteome_mediation_result <-
  nasal_microbiome_proteome_mediation_result %>%
  dplyr::mutate(treat_class = "Nasal microbiome",
                mediator_class = "Proteome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(proteome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
                phenotype_true_name = phenotype)

####lipidome
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_lipidome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_lipidome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(nasal_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

nasal_microbiome_lipidome_mediation_result <-
  mediation_result


nasal_microbiome_lipidome_mediation_result <-
  nasal_microbiome_lipidome_mediation_result %>%
  dplyr::mutate(treat_class = "Nasal microbiome",
                mediator_class = "Lipidome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(lipidome_variable_info, by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name)  %>%
dplyr::mutate(treat_true_name = Genus,
              phenotype_true_name = phenotype)


####metabolome
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_metabolome/mediation_result")
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)
load("data_analysis/mediation_analysis/sample_wise_IR/nasal_microbiome_vs_metabolome/mediation_result_inverse")

mediation_result <-
  mediation_result %>%
  dplyr::left_join(nasal_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

nasal_microbiome_metabolome_mediation_result <-
  mediation_result

nasal_microbiome_metabolome_mediation_result <-
  nasal_microbiome_metabolome_mediation_result %>%
  dplyr::mutate(treat_class = "Nasal microbiome",
                mediator_class = "Metabolome",
                phenotype_class = "Phenotype") %>%
  dplyr::left_join(metabolome_variable_info,
  by = c("mediator" = "variable_id")) %>%
  dplyr::rename(mediator_true_name = true_name) %>%
  dplyr::mutate(treat_true_name = Genus,
              phenotype_true_name = phenotype)

masstools::setwd_project()
dir.create("data_analysis/mediation_analysis/sample_wise_IR/summary", recursive = TRUE)
setwd("data_analysis/mediation_analysis/sample_wise_IR/summary")
mediation_result <-
  rbind(stool_microbiome_cytokine_mediation_result,
        stool_microbiome_proteome_mediation_result,
        stool_microbiome_lipidome_mediation_result,
        stool_microbiome_metabolome_mediation_result,
        skin_microbiome_cytokine_mediation_result,
        skin_microbiome_proteome_mediation_result,
        skin_microbiome_lipidome_mediation_result,
        skin_microbiome_metabolome_mediation_result,
        oral_microbiome_cytokine_mediation_result,
        oral_microbiome_proteome_mediation_result,
        oral_microbiome_lipidome_mediation_result,
        oral_microbiome_metabolome_mediation_result,
        nasal_microbiome_cytokine_mediation_result,
        nasal_microbiome_proteome_mediation_result,
        nasal_microbiome_lipidome_mediation_result,
        nasal_microbiome_metabolome_mediation_result
        )

dim(mediation_result)
sum(mediation_result$acme_p < 0.05)
sum(mediation_result$acme_p_inverse > 0.05)

mediation_result <- 
  mediation_result %>% 
  dplyr::filter(acme_p_inverse > 0.05 & prop_mediate < 1)
dim(mediation_result)

write.csv(mediation_result, "mediation_result_all.csv", row.names = FALSE)

#### network to show how exposome affect phenotype via omics
edge_data =
  rbind(mediation_result[,c("treat", 
                            "mediator", 
                            "treat_class",
                            "mediator_class",
                            "Genus",
                            "mediator_true_name",
                            "treat_mediator_cor",
                            "treat_mediator_p")] %>%
          dplyr::rename(from = treat, 
                        to = mediator, 
                        from_class = treat_class,
                        to_class = mediator_class,
                        from_true_name = Genus,
                        to_true_name = mediator_true_name,
                        cor = treat_mediator_cor, 
                        p = treat_mediator_p),
        mediation_result[,c("mediator", 
                            "phenotype", 
                            "mediator_class",
                            "phenotype_class",
                            "mediator_true_name",
                            "phenotype_true_name",
                            "mediator_phenotype_cor",
                            "mediator_phenotype_p")] %>%
          dplyr::rename(from = mediator, 
                        to = phenotype, 
                        from_class = mediator_class,
                        to_class = phenotype_class,
                        from_true_name = mediator_true_name,
                        to_true_name = phenotype_true_name,
                        cor = mediator_phenotype_cor, 
                        p = mediator_phenotype_p)
        # mediation_result[,c("treat", 
        #                     "phenotype", 
        #                     "treat_class",
        #                     "phenotype_class",
        #                     "treat_true_name",
        #                     "phenotype_true_name",
        #                     "treat_phenotype_cor",
        #                     "treat_phenotype_p")] %>%
        #   dplyr::rename(from = treat, 
        #                 to = phenotype, 
        #                 from_class = treat_class,
        #                 to_class = phenotype_class,
        #                 from_true_name = treat_true_name,
        #                 to_true_name = phenotype_true_name,
        #                 cor = treat_phenotype_cor, 
        #                 p = treat_phenotype_p)
  )

edge_data <- 
edge_data %>% 
  dplyr::mutate(direction = case_when(
    cor > 0 ~ "pos",
    cor < 0 ~ "neg"
  ))

node_data <-
  rbind(
    edge_data[, c("from", "from_class", "from_true_name")] %>%
      dplyr::rename(
        node = from,
        class = from_class,
        true_name = from_true_name
      ),
    edge_data[, c("to", "to_class", "to_true_name")] %>%
      dplyr::rename(
        node = to,
        class = to_class,
        true_name = to_true_name
      )
  ) %>% 
  dplyr::distinct(node, .keep_all = TRUE) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    stringr::str_detect(class, "microbiome") ~ "Treat",
                    stringr::str_detect(class, "Phenotype") ~ "Phenotype",
                    TRUE ~ "Mediator"
                  ))
  
library(tidygraph)

temp_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

library(ggraph)

my_layout <- create_layout(temp_data,
                           layout = 'linear')

my_layout$y[my_layout$class2 == "Mediator"] <- 5
my_layout$y[my_layout$class2 == "Phenotype"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

my_layout1$y[my_layout1$class2 == "Treat"] <-
  seq(from = 10, to = 90, length.out = sum(my_layout1$class2 == "Treat"))

my_layout1$y[my_layout1$class2 == "Mediator"] <-
  seq(from = 1, to = 100, length.out = sum(my_layout1$class2 == "Mediator"))

my_layout1$y[my_layout1$class2 == "Phenotype"] <-
  my_layout1$y[my_layout1$class2 == "Phenotype"] <-
  seq(from = 10, to = 90, length.out = sum(my_layout1$class2 == "Phenotype"))

plot <-
  ggraph(my_layout1) +
  geom_edge_diagonal(aes(color = direction,
                     width = -log(p, 10)),
                 show.legend = TRUE, alpha = 1) +
  geom_node_point(aes(fill = class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(
    values = c(omics_color, phenotype_color)
  ) +
  scale_color_manual(
    values = c(omics_color, phenotype_color)
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = 0,
      size = 3,
      colour = class
    ),
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class",
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  scale_edge_color_manual(values = c("pos" = "red",
                                     "neg" = viridis::cividis(n = 2)[1])) +
  ggraph::scale_edge_width(range = c(0.3, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  )

plot

ggsave(plot, filename = "mediation_analysis_all.pdf", width = 7, height = 7)


dim(edge_data)

plot(mediation_result$prop_mediate)

###
library(ggsankey)
library(dplyr)
library(ggplot2)

temp_data  <- 
mediation_result %>%
  dplyr::select(treat_class, mediator_class, phenotype_class) %>% 
    make_long(treat_class, mediator_class, phenotype_class)

plot <- 
ggplot(
  temp_data,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_manual(values = c(omics_color, phenotype_color)) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) 
plot
# ggsave(plot, filename = "sankey_plot.pdf", width = 7, height = 7)
