##avoid source
no_function()

##load data
masstools::setwd_project()
library(tidyverse)
rm(list = ls())

setwd("data_analysis/environment/")

load("expression_data")
load("sample_info")
load("variable_info")

temp_data <- 
expression_data

colnames(temp_data) <- 
  sample_info$STARTING_DATE %>% 
  stringr::str_replace("2016-", "")

rownames(temp_data) <- variable_info$true_name

temp_data <- 
temp_data %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, 
                      names_to = "date", values_to = "value") 

temp_data$variable_id <-
  factor(temp_data$variable_id, levels = c(
    "Temperature", "Humidity", "Atmospheric pressure", "Wind speed",
    "SO2", "NO2", "O3", "CO",
    "Air Quality Index", "TPM"
  ))

plot <- 
temp_data %>% 
  ggplot(aes(date, value)) +
  # geom_rect() +
  geom_line(aes(x = date, y = value, group = 1)) +
  geom_point(shape = 21,
             size = 5,
             aes(fill = variable_id), show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Temperature" = ggsci::pal_d3()(n = 10)[1], 
    "Humidity" = ggsci::pal_d3()(n = 10)[1], 
    "Atmospheric pressure" = ggsci::pal_d3()(n = 10)[1], 
    "Wind speed" = ggsci::pal_d3()(n = 10)[1], 
    "SO2" = ggsci::pal_d3()(n = 10)[3], 
    "NO2"= ggsci::pal_d3()(n = 10)[3], 
    "O3"= ggsci::pal_d3()(n = 10)[3], 
    "CO"= ggsci::pal_d3()(n = 10)[3], 
    "Air Quality Index"= ggsci::pal_d3()(n = 10)[3], 
    "TPM"= ggsci::pal_d3()(n = 10)[3]
  )) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      size = 12,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  facet_wrap(facets = vars(variable_id), scales = "free_y", nrow = 3) 
  
plot
# ggsave(plot, filename = "all_plot.pdf", width = 16, height = 8)

# for (i in 1:nrow(variable_info)) {
#   cat(i, " ")
#   temp_data <-
#     data.frame(value = as.numeric(expression_data[variable_info$variable_id[i], ]),
#                sample_info,
#                stringsAsFactors = FALSE)
#   
#   plot <-
#     temp_data %>%
#     dplyr::mutate(STARTING_DATE = stringr::str_replace(STARTING_DATE, "2016-", "")) %>%
#     ggplot(aes(x = STARTING_DATE, y = value)) +
#     geom_point(shape = 21,
#                size = 6,
#                fill = ggsci::pal_d3()(n = 10)[1]) +
#     ggrepel::geom_text_repel(aes(label = location)) +
#     theme_bw() +
#     labs(x = "", y = variable_info$variable_id[i]) +
#     theme(
#       panel.grid.minor = element_blank(),
#       axis.text.x = element_text(
#         angle = 45,
#         size = 12,
#         hjust = 1,
#         vjust = 1
#       ),
#       axis.text.y = element_text(size = 12),
#       axis.title = element_text(size = 13)
#     )
#   
#   name1 <- paste(variable_info$variable_id[i], "_1.pdf", sep = "")
#   name2 <- paste(variable_info$variable_id[i], "_2.pdf", sep = "")
#   ggsave(plot,
#          filename = name1,
#          width = 7,
#          height = 7)
#   ggsave(plot,
#          filename = name2,
#          width = 14,
#          height = 7)
# }



###PCA analysis
temp_expression_data <- 
  expression_data

temp_expression_data <- 
temp_expression_data %>% 
  apply(1, function(x){
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x
  }) %>% 
  t()

colnames(temp_expression_data) <-
  as.character(sample_info$STARTING_DATE)

#PCA analysis
###PCA for date
pca_object <- 
  prcomp(x = t(temp_expression_data), center = TRUE, scale. = TRUE)

library(wesanderson)
names(wes_palettes)
wes_palette(name = "Zissou1", n = 100, type = "continuous")

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

as.Date_origin <- function(x) {
  as.Date(x, origin = '1970-01-01')
}

plot <-
  pca_object$x %>%
  as.data.frame() %>%
  rownames_to_column(var = "date") %>%
  left_join(sample_info %>% dplyr::mutate(date = as.character(STARTING_DATE)), 
            by = c("date")) %>%
  mutate(date2 = as.character(date)) %>%
  mutate(date = as.integer(as.Date(date))) %>%
  ggplot(aes(PC1, PC2, colour = date)) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  geom_point(shape = 21,
             size = 4,
             aes(fill = date),
             color = "black") +
  guides(colour = guide_colourbar(title = "Date")) +
  scale_colour_gradientn(colours = pal,
                         labels = as.Date_origin) +
  scale_fill_gradientn(colours = pal,
                       labels = as.Date_origin) +
  ggrepel::geom_text_repel(aes(PC1, PC2, label = paste(date2, location, sep = "_")), 
                            show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent")
  ) +
  labs(
    x = paste("PC1 (", round(summary(pca_object)$importance[2, 1], 4) * 100, "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object)$importance[2, 2], 4) *
                100, "%)", sep = "")
  )

plot

# ggsave(plot,
#        file = "pca_date.pdf",
#        width = 7,
#        height = 7)

