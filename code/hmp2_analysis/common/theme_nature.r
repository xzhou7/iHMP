
library(cowplot)
library(ggplot2)
theme_nature <- function() list(
    theme_cowplot(),
    theme(
        text               = element_text(size=6),
        axis.text          = element_text(size=5),
        axis.title.x       = element_text(margin=margin(1, 0, 0.5, 0)),
        axis.title.x.top   = element_text(margin=margin(0, 0, 2, 0)),
        axis.title.y       = element_text(margin=margin(0, 1, 0, 0.5)),
        axis.title.y.right = element_text(margin=margin(0, 0, 0, 2)),
        axis.text.x        = element_text(margin=margin(1, 0, 0, 0)),
        axis.text.x.top    = element_text(margin=margin(0, 0, 1, 0)),
        axis.text.y        = element_text(margin=margin(0, 1, 0, 0)),
        axis.text.y.right  = element_text(margin=margin(0, 0, 0, 1)),
        axis.ticks         = element_line(size=0.3),
        axis.ticks.length  = unit(2, "pt"),
        axis.line          = element_line(size=0.3),
        axis.line.x        = element_line(size=0.3),
        axis.line.y        = element_line(size=0.3),
        line               = element_line(size=0.3),
        legend.margin      = margin(4, 4, 4, 4),
        legend.key.size    = unit(8, "pt"),
        legend.box.spacing = unit(4, "pt"),
        panel.spacing      = unit(1.5, "pt"),
        plot.title         = element_text(size=8),
        plot.margin        = margin(1, 1, 1, 1),
        strip.background   = element_blank(),
        strip.text         = element_text(size=6),
        strip.text.x       = element_text(margin=margin(3, 0, 3, 0)),
        strip.text.y       = element_text(margin=margin(0, 3, 0, 3))
    )
)
