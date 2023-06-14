
# Base colors

hmp2_disease_colors <- c(
	CD = "brown",
	UC = "orange",
	nonIBD = "cornflowerblue"
)

hmp2_active_disease_colors <- c(
    nonIBD = unname(hmp2_disease_colors["nonIBD"]),
    "Active CD"   = rgb(colorRamp(c(hmp2_disease_colors["CD"], "black"))(0.3), maxColorValue=255),
    "Inactive CD" = rgb(colorRamp(c(hmp2_disease_colors["CD"], "white"))(0.6), maxColorValue=255),
    "Active UC"   = rgb(colorRamp(c(hmp2_disease_colors["UC"], "black"))(0.3), maxColorValue=255),
    "Inactive UC" = rgb(colorRamp(c(hmp2_disease_colors["UC"], "white"))(0.6), maxColorValue=255)
)

# Fills/outlines

hmp2_active_disease_colors_fill <- sapply(hmp2_active_disease_colors, function(x)
    colorRampPalette(c(x, "white"))(21)[5])
hmp2_active_disease_colors_outline <- sapply(hmp2_active_disease_colors, function(x)
    colorRampPalette(c(x, "black"))(21)[5])

hmp2_disease_colors_fill <- sapply(hmp2_disease_colors, function(x)
    colorRampPalette(c(x, "white"))(21)[5])
hmp2_disease_colors_outline <- sapply(hmp2_disease_colors, function(x)
    colorRampPalette(c(x, "black"))(21)[5])
