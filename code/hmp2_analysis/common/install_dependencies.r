
source("https://bioconductor.org/biocLite.R")

# Installs packages used by the common and/or analysis scripts

install.packages(c(
	# ggplot
	"ggplot2", "reshape2", "gtable", "grid", "gridExtra",
	"ggrepel", "egg", "cowplot", "ggsignif", "ggpubr", "GGally",

	# Colors
	"viridis", "RColorBrewer",

	# Heatmaps
	"pheatmap",

	# Programming utility
	"plyr", "dplyr", "tidyr", "signal", "magrittr",
	"tibble", "stringr", "pkgmaker",

	# IO
	"data.table",

	# Ordination
	"vegan", "labdsv", "tsne", "diffusionMap", "pcaL1",

	# Modelling
	"MASS", "nlme", "lme4", "car", "metap", "ade4",

	# Optimization
	"neldermead", "seriation", "TSP",

	# Misc
	"wavethresh", "RbioRXN"
))

# Bioconductor packages
biocLite(c(
	"qvalue"
))
