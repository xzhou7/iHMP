no_source()

# set work directory
masstools::setwd_project()
library(tidyverse)
library(phyloseq)
rm(list = ls())

source("code/tools.R")

###load skin_microbiome
{
  load("data_analysis/skin_microbiome/data_preparation/sample_info")
  load("data_analysis/skin_microbiome/data_preparation/expression_data")
  load("data_analysis/skin_microbiome/data_preparation/variable_info")
  
  skin_microbiome_sample_info = sample_info
  skin_microbiome_expression_data = expression_data
  skin_microbiome_variable_info = variable_info
  
  skin_microbiome_sample_info = 
    skin_microbiome_sample_info %>% 
    dplyr::mutate(SSPG = as.numeric(SSPG)) %>% 
    dplyr::mutate(iris = case_when(
      SSPG > 125 ~ "IR",
      SSPG < 125 ~ "IS"
    )) 
  
  library(lubridate)
  
  skin_microbiome_sample_info$days = 
    skin_microbiome_sample_info$Date %>% 
    lubridate::as_date() %>% 
    yday() 
  
  skin_microbiome_sample_info = 
    skin_microbiome_sample_info %>% 
    dplyr::mutate(
      days =  yday(lubridate::as_date(Date)),
      months = lubridate::month(lubridate::as_date(skin_microbiome_sample_info$Date)),
      weeks = lubridate::week(lubridate::as_date(skin_microbiome_sample_info$Date))
    )   

  zero_percent = 
    skin_microbiome_expression_data %>% 
    apply(1, function(x){
      sum(x == 0)/ncol(skin_microbiome_expression_data)
    })
  
  sum(zero_percent > 0.99)/nrow(skin_microbiome_variable_info)
  
  ##here we remove the genus with 0 > 99%
  remain_idx = which(zero_percent < 0.99)
  
  skin_microbiome_expression_data = skin_microbiome_expression_data[remain_idx,]
  skin_microbiome_variable_info = skin_microbiome_variable_info[remain_idx,]
  
  ###only remain mike's samples
  skin_microbiome_sample_info =
    skin_microbiome_sample_info %>%
    dplyr::filter(subject_id == "69-001")
  
  skin_microbiome_expression_data = 
    skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
  
  }

##load exposome_biological
###exposome_biological
{
  load("data_analysis/exposome_biological/dna_expression_data")
  exposome_biological_expression_data <- dna_expression_data
  load("data_analysis/exposome_biological/dna_sample_info")
  exposome_biological_sample_info <- dna_sample_info
  load("data_analysis/exposome_biological/dna_variable_info")
  exposome_biological_variable_info <- dna_variable_info  
}


#######work directory
masstools::setwd_project()
setwd("data_analysis/skin_microbiome/skin_microbiome_exposome_biological")

#####match data
{
  sort(skin_microbiome_sample_info$Date)
  sort(exposome_biological_sample_info$date.start)
  
  temp1 <- 
    exposome_biological_sample_info$date.start %>% 
    as.Date() %>% 
    data.frame(date = ., 
               exposome_biological = 1,
               stringsAsFactors = FALSE) %>% 
    dplyr::arrange(date)
  
  temp2 <- 
    skin_microbiome_sample_info$Date %>% 
    as.Date() %>% 
    data.frame(date = ., 
               skin_microbiome = 1,
               stringsAsFactors = FALSE) %>% 
    dplyr::arrange(date)
  
  temp <- 
    temp1 %>% 
    dplyr::full_join(temp2, by = "date") %>% 
    arrange(date)
  
  temp <- 
    temp %>% 
    tidyr::pivot_longer(cols = -date, names_to = "class", values_to = "value")
  
  temp <- temp %>% 
    dplyr::filter(!is.na(value))
  
  diff <- 
    lapply(skin_microbiome_sample_info$Date %>% as.Date(), function(x){
      x <- abs(x - as.Date(exposome_biological_sample_info$date.start))
      as.Date(exposome_biological_sample_info$date.start)[which(x >=0 & x <= 7)]
    }) 
  
  names(diff) <- skin_microbiome_sample_info$Date %>% as.character()
  
  diff =
    diff[unlist(lapply(diff, length)) != 0]
  
  diff = diff[sort(names(diff))]
  
  ####remove some duplicated samples
  diff = diff[-1]
  
  diff$`2016-01-15` <- diff$`2016-01-15`[1]
  diff$`2016-01-19` <- diff$`2016-01-19`[2]
  diff$`2016-02-24` <- diff$`2016-02-24`[2]
  
  
  diff <- 
    purrr::map2(diff, names(diff), .f = function(x,y){
      if(length(x) == 0){
        return(NULL)
      }
      z <- data.frame(y, x, stringsAsFactors = FALSE)
      colnames(z) <- c("skin_microbiome", "exposome_biological")
      z
    }) 
  
  diff <- purrr::map2(.x = diff, .y = 1:length(diff), .f = function(x,y){
    data.frame(x, group = y, stringsAsFactors = FALSE)
  }) %>% 
    do.call(rbind, .)
  
  rownames(diff) <- NULL
  
  diff <- 
    diff %>% 
    mutate(exposome_biological = as.character(exposome_biological), 
           skin_microbiome = as.character(skin_microbiome)) %>% 
    tidyr::pivot_longer(cols = -group,
                        names_to = "class", 
                        values_to = "date")
  
  temp <-
    temp %>% 
    mutate(date = as.character(date)) %>% 
    left_join(diff, by = c("date","class"))
  
  temp$group[is.na(temp$group)] <- "No"
  
  diff
  
  library(plyr)
  
  diff <-
    diff %>% 
    plyr::dlply(.variables = .(group)) %>%
    lapply(function(x){
      x %>% 
        plyr::dlply(.variables = .(class)) %>% 
        do.call(cbind, .)
    }) %>% 
    do.call(rbind, .)
  
  diff$exposome_biological.date <- as.Date(diff$exposome_biological.date)
  
  diff$skin_microbiome.date <- as.Date(diff$skin_microbiome.date)
  
  diff$exposome_biological.group <- as.character(diff$exposome_biological.group)
  
  temp$class[temp$class == "exposome_biological"] = "Exposome"
  temp$class[temp$class == "skin_microbiome"] = "Skin microbiome"
  diff$exposome_biological.class[diff$exposome_biological.class == "exposome_biological"] = "Exposome"
  diff$skin_microbiome.class[diff$skin_microbiome.class == "skin_microbiome"] = "Skin microbiome"

}

temp_data = 
  temp %>%
  mutate(date = as.Date(date)) %>%
  dplyr::filter(date >= "2016-01-12" & date < "2016-03-10") %>%
  dplyr::mutate(class = factor(class,
                               levels = c("Exposome", "Skin microbiome"))) 

plot <-
  ggplot(aes(date, class), data = temp_data) +
  geom_segment(
    aes(
      x = exposome_biological.date,
      xend = skin_microbiome.date,
      y = exposome_biological.class,
      yend = skin_microbiome.class
    ),
    color = "black",
    data = diff,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = class),
    color = "black",
    alpha = 1,
    size = 4,
    shape = 21,
    show.legend = FALSE,
    data = temp_data
  ) +
  scale_fill_manual(values = omics_color) +
  scale_x_continuous(
    trans = "date",
    breaks = c(as.Date(temp$date)),
    labels = as.character(temp$date)
  ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "exposome_biological_skin_microbiome_match.pdf", width = 10, height = 7)


# #-------------------------------------------------------------------------------
# ###prepare data
# dim(skin_microbiome_expression_data)
# dim(skin_microbiome_sample_info)
# 
# skin_microbiome_sample_info <-
#   skin_microbiome_sample_info %>%
#   dplyr::filter(Date %in% unique(diff$skin_microbiome.date))
# 
# skin_microbiome_sample_info
# 
# skin_microbiome_expression_data <-
#   skin_microbiome_expression_data %>%
#   dplyr::select(one_of(skin_microbiome_sample_info$sample_id))
# 
# skin_microbiome_sample_info <-
#   skin_microbiome_sample_info %>%
#   dplyr::left_join(diff[,c("skin_microbiome.date", "skin_microbiome.group")] %>% dplyr::distinct(skin_microbiome.date, skin_microbiome.group),
#                    by = c("Date" = "skin_microbiome.date")) %>%
#   dplyr::rename(group = skin_microbiome.group)
# 
# exposome_biological_sample_info <-
#   exposome_biological_sample_info %>%
#   dplyr::filter(as.Date(date.start) %in% diff$exposome_biological.date)
# 
# exposome_biological_expression_data <-
#   exposome_biological_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "Sample_ID") %>%
#   dplyr::mutate(Sample_ID = stringr::str_replace(Sample_ID, "_[0-9]{1,2}", "")) %>%
#   plyr::dlply(.variables = "Sample_ID") %>%
#   lapply(function(x){
#     apply(x[,-1], 2, mean)
#   }) %>%
#   do.call(rbind, .) %>%
#   t() %>%
#   as.data.frame()
# 
# ###browser
# exposome_biological_expression_data <-
#   exposome_biological_expression_data %>%
#   dplyr::select(one_of(exposome_biological_sample_info$sample_id))
# 
# exposome_biological_sample_info <-
#   exposome_biological_sample_info %>%
#   mutate(date.start = as.Date(date.start)) %>%
#   dplyr::left_join(diff[,c("exposome_biological.date", "exposome_biological.group")], by = c("date.start" = "exposome_biological.date")) %>%
#   dplyr::rename(group = exposome_biological.group)
# 
# dim(skin_microbiome_sample_info)
# dim(exposome_biological_sample_info)
# 
# ###sort date
# skin_microbiome_sample_info$sample_id == colnames(skin_microbiome_expression_data)
# exposome_biological_sample_info$sample_id == colnames(exposome_biological_expression_data)
# 
# skin_microbiome_sample_info =
# skin_microbiome_sample_info %>%
#   dplyr::arrange(Date)
# 
# skin_microbiome_expression_data =
#   skin_microbiome_expression_data[,skin_microbiome_sample_info$sample_id]
# 
# exposome_biological_sample_info =
#   exposome_biological_sample_info %>%
#   dplyr::arrange(date.start)
# 
# exposome_biological_expression_data =
#   exposome_biological_expression_data[,exposome_biological_sample_info$sample_id]
# 
# exposome_biological_sample_info =
# exposome_biological_sample_info %>%
#   dplyr::mutate(date = date.start)
# 
# skin_microbiome_sample_info =
#   skin_microbiome_sample_info %>%
#   dplyr::mutate(date = Date)
# 
# save(exposome_biological_expression_data, file = "exposome_biological_expression_data")
# save(skin_microbiome_expression_data, file = "skin_microbiome_expression_data")
# 
# save(exposome_biological_sample_info, file = "exposome_biological_sample_info")
# save(skin_microbiome_sample_info, file = "skin_microbiome_sample_info")



####correlation analysis
load("exposome_biological_expression_data")
load("skin_microbiome_expression_data")
load("exposome_biological_sample_info")
load("skin_microbiome_sample_info")

##remove microbiome with a lot of zero values
zero_percent = 
  skin_microbiome_expression_data %>% 
  apply(1, function(x){
    sum(x == 0)/ncol(skin_microbiome_expression_data)
  })

plot(zero_percent)

sum(zero_percent < 0.5)

remain_idx = which(zero_percent < 0.5)

length(remain_idx)

skin_microbiome_expression_data = 
  skin_microbiome_expression_data[remain_idx,]

skin_microbiome_variable_info =
  skin_microbiome_variable_info[remain_idx,]


##remove microbiome with a lot of zero values
zero_percent = 
  exposome_biological_expression_data %>% 
  apply(1, function(x){
    sum(x == 0)/ncol(exposome_biological_expression_data)
  })

plot(zero_percent)

sum(zero_percent < 0.5)

remain_idx = which(zero_percent < 0.5)

length(remain_idx)

exposome_biological_expression_data = 
  exposome_biological_expression_data[remain_idx,]

exposome_biological_variable_info =
  exposome_biological_variable_info[remain_idx,]

skin_microbiome_sample_info$sample_id == colnames(skin_microbiome_expression_data)
colnames(skin_microbiome_expression_data) <- as.character(skin_microbiome_sample_info$date)
colnames(exposome_biological_expression_data) <- colnames(skin_microbiome_expression_data)

# #######correlation analysis
# exposome_biological_expression_data <- log(exposome_biological_expression_data + 1, 2)
# 
# dim(exposome_biological_expression_data)
# dim(skin_microbiome_expression_data)
# 
# ###correct fiber for skin_microbiome
# skin_microbiome_expression_data1 <-
#   purrr::map(as.data.frame(t(skin_microbiome_expression_data)), .f = function(x){
#     temp_data <-
#       data.frame(fiber = c(1,
#                            1,
#                            1,
#                            0,
#                            1,
#                            1),
#                  x,
#                  stringsAsFactors = FALSE)
# 
#     lm_result <- lm(formula = x ~ fiber, data = temp_data)
#     lm_result$residuals
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# temp_data <-
#   skin_microbiome_expression_data1
# 
# colnames(temp_data) <-
#   colnames(skin_microbiome_expression_data1) <-
#   colnames(skin_microbiome_expression_data)
# 
# 
# ###calculate correlation between skin_microbiome and exposome
# exposome_biological_expression_data1 <- exposome_biological_expression_data
# save(exposome_biological_expression_data1, file = "exposome_biological_expression_data1")
# save(skin_microbiome_expression_data1, file = "skin_microbiome_expression_data1")
# 
# cor_value <-
#   cor(x = t(as.matrix(exposome_biological_expression_data)),
#       y = t(as.matrix(skin_microbiome_expression_data1)),
#       method = "spearman")
# 
# cor_value <-
#   cor_value %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "from") %>%
#   tidyr::pivot_longer(-from, names_to = "to", values_to = "cor")
# 
# library(plyr)
# 
# p_value <-
#   purrr::map(as.data.frame(t(cor_value)), .f = function(x){
#     value1 <- as.numeric(exposome_biological_expression_data[x[1],])
#     value2 <- as.numeric(skin_microbiome_expression_data1[x[2],])
#     cor.test(value1, value2, method = "spearman")$p.value
#   }) %>%
#   unlist()
# 
# cor_value <-
#   data.frame(cor_value, p_value, stringsAsFactors = FALSE)
# 
# plot(density(cor_value$p_value))
# 
# library(plyr)
# cor_value <-
#   cor_value %>%
#   plyr::dlply(.variables = .(from)) %>%
#   purrr::map(.f = function(x){
#     # x <- x %>%
#       # dplyr::filter(abs(cor) > 0.9)
#     fdr <- p.adjust(x$p_value, method = "fdr")
#     x <-
#       data.frame(x, fdr, stringsAsFactors = FALSE)
#     x
#   })
# 
# cor_value <-
#   cor_value %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# # cor_value <-
# #   cor_value %>%
# #   dplyr::filter(abs(cor) > 0.9 & cor != 1)
# 
# dim(cor_value)
# 
# save(cor_value, file = "cor_value")
load('cor_value')

cor_value1 <- 
  cor_value %>% 
  dplyr::filter(fdr < 0.05)

dim(cor_value1)

cor_value1 <-
  cor_value1 %>%
  dplyr::left_join(exposome_biological_variable_info[,c(1:2)],
                   by = c("from" = "variable_id")) %>%
  dplyr::mutate(exposome_biological_id = from) %>%
  dplyr::left_join(skin_microbiome_variable_info, by = c("to" = "variable_id")) %>%
  dplyr::mutate(skin_microbiome_id = to) %>%
  dplyr::filter(!is.na(skin_microbiome_id))

dir.create("scatter_plot")

for(idx in 1:nrow(cor_value1)) {
  cat(idx, " ")
  path1 <- file.path(cor_value1$from[idx])
  temp_data <-
    data.frame(
      date = as.character(skin_microbiome_sample_info$date),
      exp = as.numeric(exposome_biological_expression_data[cor_value1$from[idx], ]),
      pro = as.numeric(skin_microbiome_expression_data1[cor_value1$to[idx], ]),
      stringsAsFactors = FALSE
    )
  plot <-
    temp_data %>%
    ggplot(aes(exp, pro)) +
    geom_point() +
    geom_smooth(method = "lm", color = "skyblue") +
    ggrepel::geom_label_repel(aes(x = exp, pro, label = date)) +
    labs(
      x = paste("Exposome biological: ", cor_value1$exposome_biological_id[idx], sep = ""),
      y = paste("Skin microbiome: " , cor_value1$skin_microbiome_id[idx]),
      sep = ""
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    ) +
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      label = paste(
        "Correlation: ",
        round(cor_value1$cor[idx], 2),
        "\nFDR adjusted P value: ",
        round(cor_value1$fdr[idx], 3),
        sep = ""
      ),
      vjust = 2,
      hjust = -1
    )
  
  name <- paste(cor_value1$from[idx], "_",
                cor_value1$to[idx], ".pdf", sep = "")
  
  dir.create(file.path("scatter_plot",path1))
  
  ggsave(
    plot,
    filename = file.path("scatter_plot",path1, name),
    width = 7,
    height = 7,
    bg = "transparent"
  )
  
}

cor_value1$from %>% unique()

###correlation network for pro and exp
cor_value1$from %>% unique() %>% length()
cor_value1$to %>% unique() %>% length()

library(igraph)
library(ggraph)
library(tidygraph)

###network for all the exposome and skin_microbiome
edge_data <-  
  cor_value1 %>% 
  # dplyr::filter(from %in% cluster1) %>%
  dplyr::rename(from = from, 
                to = to, 
                Correlation = cor) %>% 
  dplyr::mutate(fdr = -log(fdr, 10))

node_data <- 
  cor_value1 %>% 
  # dplyr::filter(from %in% cluster1) %>%
  dplyr::rename(from = from, to = to) %>% 
  dplyr::select(from, to) %>% 
  tidyr::pivot_longer(cols = c(from, to), 
                      names_to = "class", values_to = "node") %>% 
  dplyr::mutate(class1 = case_when(
    stringr::str_detect(class, "from") ~ "Exposome",
    TRUE ~ "Skin microbiome"
  )) %>% 
  dplyr::select(node, class1) %>% 
  dplyr::rename(Class = class1) %>%
  dplyr::distinct(node, .keep_all = TRUE)

node_data <- 
  node_data %>% 
  dplyr::arrange(Class)

# node_data <- 
#   node_data %>% 
#   dplyr::left_join(exposome_biological_variable_info[,c("peak_ID", "MetabID")], by = c("node" = "peak_ID")) %>% 
#   dplyr::mutate(compound.name = case_when(
#     !is.na(MetabID) ~ MetabID,
#     TRUE ~ node
#   )) %>% 
#   dplyr::select(node, Class, compound.name)

node_data %>% 
  dplyr::left_join(skin_microbiome_variable_info)

temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = TRUE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

plot1 <-
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(aes(color = Correlation),
                show.legend = TRUE) +
  geom_node_point(aes(fill = Class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = omics_color) +
  scale_color_manual(values = omics_color) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = node,
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      size = 3,
      color = Class
    ),
    size = 3,
    alpha = 1, 
    show.legend = FALSE
  ) +
  guides(edge_width = guide_legend(title = "-log10(FDR adjusted P value)", 
                                   override.aes = list(shape = NA)),
         edge_color = ggraph::guide_edge_colorbar(title = "Spearman correlation"),
         fill = guide_legend(title = "Class", 
                             override.aes = list(size = 4, linetype = "blank")),
         size = guide_legend(title = "Degree", override.aes = list(linetype = 0))) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.2, 2)) +
  scale_size_continuous(range = c(1, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

plot1

ggsave(
  plot1,
  filename = "exposome_biological_skin_microbiome_correlation_network.pdf",
  width = 8.5,
  height = 7,
  bg = "transparent"
)
