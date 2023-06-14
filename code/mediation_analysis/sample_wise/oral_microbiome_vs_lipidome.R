####reference 
####https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171

no_function()
# set work directory

setwd(masstools::get_project_wd())
library(tidyverse)
rm(list = ls())

source("code/tools.R")

{
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/oral_microbiome_expression_data"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/oral_microbiome_variable_info"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/oral_microbiome_sample_info"
  )
  
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/lipidome_expression_data"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/lipidome_variable_info"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/lipidome_sample_info"
  )
}


{
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_phenotype_sample_wise/phenotype_expression_data"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_phenotype_sample_wise/phenotype_sample_info"
  )
  load(
    "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_phenotype_sample_wise/phenotype_variable_info"
  )
}

colnames(oral_microbiome_expression_data)
colnames(lipidome_expression_data)
colnames(phenotype_expression_data)

intersect_sample_id <-
  Reduce(f = intersect, x = list(
    colnames(oral_microbiome_expression_data),
    colnames(lipidome_expression_data),
    colnames(phenotype_expression_data)
  ))

oral_microbiome_expression_data <- 
  oral_microbiome_expression_data[,intersect_sample_id]

lipidome_expression_data <- 
  lipidome_expression_data[,intersect_sample_id]

phenotype_expression_data <- 
  phenotype_expression_data[,intersect_sample_id]

oral_microbiome_sample_info <- 
  oral_microbiome_sample_info[match(intersect_sample_id, oral_microbiome_sample_info$sample_id),]

lipidome_sample_info <- 
  lipidome_sample_info[match(intersect_sample_id, lipidome_sample_info$sample_id),]

phenotype_sample_info <- 
  phenotype_sample_info[match(intersect_sample_id, phenotype_sample_info$sample_id),]

#####remove some genus
##remove Archaea
remove_name <-
  c(
    "Zea",
    oral_microbiome_variable_info %>%
      dplyr::filter(Kingdom == "Archaea") %>%
      dplyr::pull(Genus)
  ) %>%
  unique()

if (length(remove_name) > 0) {
  oral_microbiome_variable_info <-
    oral_microbiome_variable_info %>%
    dplyr::filter(!Genus %in% remove_name)
  
  oral_microbiome_expression_data <-
    oral_microbiome_expression_data[oral_microbiome_variable_info$variable_id, ]
}

dim(oral_microbiome_expression_data)
dim(lipidome_expression_data)

###finally, for stool microbiome, 107 genus, for lipidome, 302 proteins
###563 samples

######--------------------------------------------------------------------------
##for raw data, we just log(x+1, 2)
library(plyr)

###because our microbiome are percentage data, so here we use the CTL method
library(compositions)
oral_microbiome_expression_data =
  oral_microbiome_expression_data %>%
  purrr::map(function(x) {
    x = compositions::clr(x) %>%
      as.numeric()
    x
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

rownames(oral_microbiome_expression_data) = oral_microbiome_variable_info$variable_id

#####
dim(lipidome_expression_data)
dim(oral_microbiome_expression_data)

oral_microbiome_sample_info$subject_id == lipidome_sample_info$subject_id

load(
  "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_lipidome/oral_microbiome_lipidome_lm_adjusted_cor_spearman"
)

oral_microbiome_lipidome_cor <- 
  oral_microbiome_lipidome_lm_adjusted_cor_spearman %>% 
  dplyr::filter(p_adjust < 0.2)

load(
  "data_analysis/correlation_network/whole_data_set/oral_microbiome_vs_phenotype_sample_wise/oral_microbiome_phenotype_lm_adjusted_cor_spearman"
)

oral_microbiome_phenotype_cor <- 
  oral_microbiome_phenotype_lm_adjusted_cor_spearman %>% 
  dplyr::filter(p_adjust < 0.2)


load(
  "data_analysis/correlation_network/whole_data_set/lipidome_vs_phenotype_sample_wise/lipidome_phenotype_lm_adjusted_cor_spearman"
)

lipidome_phenotype_cor <- 
  lipidome_phenotype_lm_adjusted_cor_spearman %>% 
  dplyr::filter(p_adjust < 0.2)

dim(oral_microbiome_lipidome_cor)
dim(oral_microbiome_phenotype_cor)

unique(oral_microbiome_lipidome_cor$microbiome)
unique(oral_microbiome_phenotype_cor$microbiome)

intersect_id <- 
  intersect(unique(oral_microbiome_lipidome_cor$microbiome),
            unique(oral_microbiome_phenotype_cor$microbiome))

oral_microbiome_lipidome_cor <- 
  oral_microbiome_lipidome_cor %>% 
  dplyr::filter(microbiome %in% intersect_id)

oral_microbiome_phenotype_cor <- 
  oral_microbiome_phenotype_cor %>% 
  dplyr::filter(microbiome %in% intersect_id)

######work directory
setwd(masstools::get_project_wd())
dir.create("data_analysis/mediation_analysis/sample_wise/oral_microbiome_vs_lipidome",
           recursive = TRUE)
setwd("data_analysis/mediation_analysis/sample_wise/oral_microbiome_vs_lipidome")

###here we use the lm_adjusted_cor
sum(oral_microbiome_lipidome_cor$p_adjust < 0.2)

library(mediation)

plot(density(as.numeric(lipidome_expression_data[,1])))

mediation_result = NULL

for(i in 1:nrow(oral_microbiome_lipidome_cor)) {
  cat(i, "\n")
  microbiome_id = oral_microbiome_lipidome_cor$microbiome[i]
  internal_ome_id = oral_microbiome_lipidome_cor$metabolite[i]

  temp_data <-
  data.frame(microbiome = as.numeric(oral_microbiome_expression_data[microbiome_id,]),
             internal_ome = as.numeric(lipidome_expression_data[internal_ome_id,]),
             oral_microbiome_sample_info
             )

  temp_data <-
  temp_data %>%
    dplyr::mutate(Gender = case_when(
      Gender == "F" ~ 0,
      Gender == "M" ~ 1
    ))

  glm_reg1 =
    glm(
      internal_ome ~ microbiome + Gender + Adj.age + BMI,
      family = gaussian,
      temp_data
    )

  p <- broom::tidy(glm_reg1) %>%
    as.data.frame()

  if(p$p.value[2] > 0.05){
    next()
  }

  temp_oral_microbiome_phenotype_cor =
    oral_microbiome_phenotype_cor %>%
    dplyr::filter(microbiome %in% microbiome_id)

  if(nrow(temp_oral_microbiome_phenotype_cor) == 0){
    next()
  }

  for (j in 1:nrow(temp_oral_microbiome_phenotype_cor)) {
    cat(j, " ")
    phenotype_id <-
      temp_oral_microbiome_phenotype_cor$metabolite[j]

    temp_data2 =
      data.frame(phenotype = as.numeric(phenotype_expression_data[phenotype_id,]),
                 temp_data)

    temp_data2$phenotype[is.na(temp_data2$phenotype)] <- 0

    glm_reg2 =
      glm(
        phenotype ~ microbiome + internal_ome + Gender + Adj.age + BMI,
        family = gaussian,
        temp_data2
      )

    result = mediate(
      model.m = glm_reg1,
      model.y = glm_reg2,
      treat = "microbiome",
      mediator = "internal_ome",
      boot = TRUE,
      sims = 500
    )

    return_result <- data.frame(
      phenotype = phenotype_id,
      mediator = internal_ome_id,
      treat = microbiome_id,
      acme = result$d1,
      acme_ci_lower = result$d1.ci[1],
      acme_ci_upper = result$d1.ci[2],
      acme_p = result$d1.p,
      ade = result$z1,
      ade_ci_lower = result$z1.ci[1],
      ade_ci_upper = result$z1.ci[2],
      ade_p = result$z1.p,
      total_effect = result$tau.coef,
      total_effect_ci_lower = result$tau.ci[1],
      total_effect_ci_upper = result$tau.ci[2],
      total_effect_p = result$tau.p,
      prop_mediate = result$n1,
      prop_mediate_ci_lower = result$n1.ci[1],
      prop_mediate_ci_upper = result$n1.ci[2],
      prop_mediate_p = result$n1.p
    )
  mediation_result = rbind(mediation_result, return_result)
  }
}

rownames(mediation_result) <- NULL

mediation_result <-
mediation_result %>%
  dplyr::left_join(oral_microbiome_lipidome_lm_adjusted_cor_spearman %>% dplyr::select(-p_adjust2),
                   by = c("treat" = "microbiome",
                          "mediator" = "metabolite")) %>%
  dplyr::rename(treat_mediator_cor = cor,
                treat_mediator_p = p,
                treat_mediator_p_adjust = p_adjust) %>%
  dplyr::left_join(oral_microbiome_phenotype_lm_adjusted_cor_spearman %>% dplyr::select(-p_adjust2),
                   by = c("treat" = "microbiome",
                          "phenotype" = "metabolite")) %>%
  dplyr::rename(treat_phenotype_cor = cor,
                treat_phenotype_p = p,
                treat_phenotype_p_adjust = p_adjust) %>%
  dplyr::left_join(lipidome_phenotype_lm_adjusted_cor_spearman %>% dplyr::select(-p_adjust2),
                   by = c("mediator" = "microbiome",
                          "phenotype" = "metabolite")) %>%
  dplyr::rename(mediator_phenotype_cor = cor,
                mediator_phenotype_p = p,
                mediator_phenotype_p_adjust = p_adjust)

save(mediation_result, file = "mediation_result")
load("mediation_result")
dim(mediation_result)
mediation_result <-
  mediation_result %>%
  dplyr::filter(prop_mediate > 0 & acme_p < 0.05)
dim(mediation_result)

###inverse mediation effect
mediation_result_inverse <- NULL
for(i in 1:nrow(mediation_result)){
  cat(i, " ")
  x <- as.character(mediation_result[i,1:3])
  microbiome_id = x[3]
  internal_ome_id = x[2]
  phenotype_id = x[1]

    temp_data <-
    data.frame(microbiome = as.numeric(oral_microbiome_expression_data[microbiome_id,]),
               internal_ome = as.numeric(lipidome_expression_data[internal_ome_id,]),
               phenotype = as.numeric(phenotype_expression_data[phenotype_id,]),
               oral_microbiome_sample_info
               ) %>%
      dplyr::filter(!is.na(phenotype))

    temp_data <-
    temp_data %>%
      dplyr::mutate(Gender = case_when(
        Gender == "F" ~ 0,
        Gender == "M" ~ 1
      ))

    glm_reg1 =
      glm(
        phenotype ~ microbiome + Gender + Adj.age + BMI,
        family = gaussian,
        temp_data
      )

    glm_reg2 =
      glm(internal_ome ~ microbiome + phenotype + Gender + Adj.age + BMI,
          family = gaussian,
          temp_data)

        result = mediate(
          model.m = glm_reg1,
          model.y = glm_reg2,
          treat = "microbiome",
          mediator = "phenotype",
          boot = TRUE,
          sims = 500
        )

        return_result <- data.frame(
          phenotype = phenotype_id,
          mediator = internal_ome_id,
          treat = microbiome_id,
          acme = result$d1,
          acme_p = result$d1.p
        ) %>%
          dplyr::rename(acme_inverse = acme,
                        acme_p_inverse = acme_p)
        mediation_result_inverse <-
          rbind(mediation_result_inverse,
                return_result)

}

save(mediation_result_inverse, file = "mediation_result_inverse")
load("mediation_result_inverse")

##add microbiome information
mediation_result <-
  mediation_result %>%
  dplyr::left_join(oral_microbiome_variable_info, by = c("treat" = "variable_id")) %>%
  dplyr::select(treat, mediator, phenotype, Kingdom:Species, everything())

mediation_result <-
  mediation_result %>%
  dplyr::left_join(mediation_result_inverse,
                   by = c("phenotype", "mediator", "treat"))

#### network to show how exposome affect phenotypes via omics
edge_data =
  rbind(
    mediation_result[, c("treat", "mediator", "treat_mediator_cor", "treat_mediator_p")] %>%
      dplyr::rename(
        from = treat,
        to = mediator,
        cor = treat_mediator_cor,
        p = treat_mediator_p
      ),
    mediation_result[, c("mediator",
                         "phenotype",
                         "mediator_phenotype_cor",
                         "mediator_phenotype_p")] %>%
      dplyr::rename(
        from = mediator,
        to = phenotype,
        cor = mediator_phenotype_cor,
        p = mediator_phenotype_p
      )
  )

node_data =
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>%
  dplyr::mutate(
    class1 = case_when(
      node %in% oral_microbiome_variable_info$variable_id ~ "Oral microbiome",
      node %in% lipidome_variable_info$variable_id ~ "Lipidome",
      node %in% phenotype_variable_info$variable_id ~ "Phenotype"
    )
  ) %>%
  dplyr::select(node, class1) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::left_join(oral_microbiome_variable_info[, c("variable_id", "Genus")],
                   by = c("node" = "variable_id")) %>%
  dplyr::mutate(true_name = case_when(!is.na(Genus) ~ Genus,
                                      TRUE ~ node)) %>%
  dplyr::select(-Genus)

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

my_layout$y[my_layout$class == "Lipidome"] <- 5
my_layout$y[my_layout$class == "Phenotype"] <- 10

my_layout1 <-
  my_layout

my_layout1$x <-
  my_layout$y

my_layout1$y <-
  my_layout$x

my_layout1$y[my_layout1$class == "Oral microbiome"] <-
  seq(
    from = 1,
    to = 100,
    length.out = sum(my_layout1$class == "Oral microbiome")
  )

my_layout1$y[my_layout1$class == "Lipidome"] <-
  seq(
    from = 1,
    to = 100,
    length.out = sum(my_layout1$class == "Lipidome")
  )

my_layout1$y[my_layout1$class == "Phenotype"] <-
  my_layout1$y[my_layout1$class == "Phenotype"] <-
  seq(
    from = 30,
    to = 70,
    length.out = sum(my_layout1$class == "Phenotype")
  )

plot <-
  ggraph(my_layout1) +
  geom_edge_link(aes(color = cor,
                     width = -log(p, 10)),
                 show.legend = TRUE,
                 alpha = 0.5) +
  geom_node_point(aes(fill = class1,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = c(omics_color[c("Oral microbiome", "Cytokine")], phenotype_color)) +
  scale_color_manual(values = c(omics_color[c("Oral microbiome", "Cytokine")], phenotype_color)) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1,
      label = true_name,
      hjust = ifelse(class1 == "Lipidome", 1, 0),
      size = 3,
      colour = class1
    ),
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
    fill = guide_legend(
      title = "Class",
      override.aes = list(size = 4, linetype = "blank")
    ),
    size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
  ) +
  scale_edge_color_gradient2(low = viridis::cividis(n = 2)[1],
                             mid = "white",
                             high = "red") +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    # legend.position = c(1,0), legend.justification = c(1,0)
  )

plot

ggsave(plot, filename = "mediation_analysis.pdf", width = 7, height = 7)


###some examples
length(unique(mediation_result$phenotype))
length(unique(mediation_result$mediator))
length(unique(mediation_result$treat))

write.csv(mediation_result, "mediation_result.csv", row.names = FALSE)
