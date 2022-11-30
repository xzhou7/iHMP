##avoid source
no_function()

masstools::setwd_project()
rm(list=ls())

library(tidyverse)

load("data_analysis/exposome_chemical/sample_info")
load("data_analysis/exposome_chemical/variable_info")
load("data_analysis/exposome_chemical/expression_data")

# sample_info <-
#   sample_info %>% 
#   dplyr::filter(start_date >= as.Date("2016-01-15") & start_date <= as.Date("2016-03-03"))

sample_id <-
  sample_info$sample_id

sample_id <-
  c(paste(sample_id, c(1), sep = "_"),
    paste(sample_id, c(2), sep = "_"),
    paste(sample_id, c(3), sep = "_"))

expression_data <-
  expression_data %>% 
  dplyr::select(one_of(sample_id))

setwd("data_analysis/exposome_chemical/")

# openxlsx::write.xlsx(sample_info, file = "exposome_sample_info.xlsx")
# openxlsx::write.xlsx(expression_data, file = "exposome_expression_data.xlsx")
# openxlsx::write.xlsx(variable_info, file = "exposome_variable_info.xlsx")

library("ggmap")

lat_lon_df <-
  sample_info %>%
  dplyr::select(location,
                latitude,
                longitude) %>%
  unique() %>%
  ungroup() %>%
  dplyr::rename(id = location) %>%
  mutate(id = factor(id))

lat_lon_df <- 
  lat_lon_df %>% 
  distinct(latitude, longitude, .keep_all = TRUE)

# Library
library(ggmap)
library(gridExtra)

# For stamen map, you have to give the extremity of the window you are looking at. here is an example with the watercolor background (Around brisbane)
us_map <-
  get_stamenmap(
    bbox = c(-127, 23, -92, 47),
    zoom = 6,
    maptype = "watercolor"
    # color = "color"
  )

plot <- 
ggmap(us_map) +
  geom_point(
    aes(x = longitude, 
        y = latitude),
    data = lat_lon_df,
    fill = ggsci::pal_aaas()(n=10)[2],
    size = 6, 
    shape = 21,
    alpha = 0.8
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(colour = "orange"),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 2
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  ggrepel::geom_text_repel(mapping = aes(x = longitude,
                                         y = latitude,
                                         label = id),
                           data = lat_lon_df, 
                           color = ggsci::pal_aaas()(n=10)[2],
                           size = 6
  )

plot

# ggsave(plot, filename = "sample_collection_us_map.pdf", width = 9, height = 7)

cal_map <- get_stamenmap(
  bbox = c(-122.7, 37.3, -122, 38.25),
  zoom = 12,
  maptype = "watercolor"
)

plot <- 
  ggmap(cal_map) +
  geom_point(
    aes(x = longitude, 
        y = latitude),
    data = lat_lon_df,
    fill = ggsci::pal_aaas()(n=10)[2],
    size = 6, 
    shape = 21,
    alpha = 0.8
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(colour = "orange"),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 2
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  ggrepel::geom_text_repel(mapping = aes(x = longitude,
                                         y = latitude,
                                         label = id),
                           data = lat_lon_df, 
                           color = ggsci::pal_aaas()(n=10)[2],
                           size = 6
  )

plot

# ggsave(plot, filename = "sample_collection_ca_map.pdf", width = 9, height = 7)

plot1 <- 
sample_info %>% 
  mutate(start_date = as.Date(start_date)) %>% 
  ggplot(aes(x = location)) +
  geom_point(aes(x = start_date, 
                 y = location, 
                 color = location),
             show.legend = FALSE,
             shape = 16,
             alpha = 0.7,
             size = 4) +
  labs(x = "", y = "") +
  ggsci::scale_color_aaas() +
  scale_x_continuous(trans = "date",
                     breaks = c(as.Date(sample_info$start_date)),
                     labels = as.character(sample_info$start_date)
  ) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, 
                                   angle = 45,
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = ggsci::pal_aaas()(length(unique(sample_info$location)))),
        # panel.grid = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        panel.grid.minor = element_blank()
  )

plot1

plot2 <-
  sample_info %>% 
  mutate(start_date = as.Date(start_date)) %>% 
  ggplot(aes(x = start_date)) +
  labs(x = "", y = "Sample number") +
  geom_histogram(fill = "#CCEBC5") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    panel.grid.minor = element_blank()
  )

plot2

plot3 <- 
  sample_info %>% 
  ggplot(aes(x = location)) +
  geom_bar(aes(fill = location), show.legend = FALSE) +
  labs(x = "", y = "Sample number") +
  ggsci::scale_fill_aaas() +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip()

plot3

library(patchwork)

plot <-
  {
    plot2 + plot2 + plot_layout(ncol = 2, widths = c(4, 1))
  } -
  {
    plot1 + plot3 + plot_layout(ncol = 2, widths = c(4, 1))
  } +
  plot_layout(ncol = 1, heights = c(1, 4))

plot

# ggsave(plot, file = "sample_collection.pdf", width = 9, height = 7)

plot3 <- 
  sample_info %>% 
  dplyr::group_by(location) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(aes(x = location, y = n)) +
  # geom_line(aes(group = 1)) +
  geom_segment(aes(color = location, 
                   x = location, xend = location, y = 0, yend = n),
               show.legend = FALSE) +
  geom_point(shape = 21, 
             aes(fill = location), show.legend = FALSE,
             size = 5) +
  labs(x = "", y = "Sample number") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = c(0,2,4,6,8,9), labels = c(0,2,4,6,8,9)) +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.text.y = element_text(angle = 90, vjust = 1, size = 12),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip()

plot3

# ggsave(plot3, filename = "location_sample.pdf", width = 2.5, height = 7)

###data qualitiy
dim(expression_data)
colnames(expression_data)
dim(sample_info)

###correct location and date
temp_data <- 
  expression_data %>% 
  dplyr::select(-contains("Blank"))

subject_id <- 
  colnames(temp_data) %>% 
  stringr::str_replace("\\_[0-9]{1}", "")

temp_sample_info <- 
data.frame(sample_id = colnames(temp_data),
           subject_id,
           stringsAsFactors = FALSE
           ) %>% 
  dplyr::left_join(sample_info, by = c("subject_id" = "sample_id"))

temp_sample_info$time <-
  temp_sample_info$start_date %>% 
  as.Date() %>% 
    `-`(as.Date("2016-01-01")) %>% 
    as.numeric()

colnames(temp_data) == temp_sample_info$sample_id

temp_data[1,]

###log
temp_data[1,] %>% 
  as.numeric() %>% 
  density() %>% 
  plot()

temp_data <- log(as.matrix(temp_data) + 1, 2)

###correct location and time
temp_data <- 
apply(temp_data, 1, function(y){
  y <- as.numeric(y)
  temp <- data.frame(y, location = temp_sample_info$location, 
                     time = temp_sample_info$time,
                     stringsAsFactors = FALSE)
  lm_res <- lm(formula = y ~ time + location, data = temp)
  lm_res$residuals
}) %>% 
  t() %>% 
  as.data.frame()

colnames(temp_data) <- temp_sample_info$sample_id

pca_object <- prcomp(t(temp_data),
                     scale. = TRUE, center = TRUE)


x <- pca_object$x[,1:2]

x <- data.frame(x, temp_sample_info, stringsAsFactors = FALSE)

library(ggfortify)
library(cowplot)
library(ggforce)
x$start_date <- as.character(x$start_date)

plot <- 
autoplot(object = pca_object,
         data = x, 
         colour = "start_date", 
         shape = "start_date",
         size = 8,
         alpha = 1,
         frame = TRUE, 
         # frame.type = 'norm',
         variance_percentage = TRUE,
         scale = 0) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  theme_bw() +
  scale_shape_manual(values = c(1:21)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0)) 
  # facet_zoom(xlim = c(-5, 5)
  #            # ylim = c(-2, 2)
  #            ) 

plot

# ggsave(plot, file = "pca_data_quality.pdf", width = 7, height = 7)











