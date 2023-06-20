#heatmap for the cutokine richness plot 

library(ComplexHeatmap)
robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

library(reshape2)
library(dplyr)
library(circlize)
library(dendextend)
library(ggplot2)
library(stringr)

setwd("~/Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/Analysis/CytokineRichness/")
stool.result <- read.csv("./stool.richness.20.csv", header = T, row.names = 1)
skin.result <- read.csv("./skin.richness.20.csv", header = T, row.names = 1)
oral.result <- read.csv("./oral.richness.20.csv", header = T, row.names = 1)
nasal.result <- read.csv("./nasal.richness.20.csv", header = T, row.names = 1)

#remove control from cytokine data (if there are such data, currently we remove it from the modelling step)
stool.result <- filter(stool.result, !str_detect(stool.result$Cytokine, "CHEX"))
skin.result <- filter(skin.result, !str_detect(skin.result$Cytokine, "CHEX"))
oral.result <- filter(oral.result, !str_detect(oral.result$Cytokine, "CHEX"))
nasal.result <- filter(nasal.result, !str_detect(nasal.result$Cytokine, "CHEX"))

#Process stool data
wide.stool <- dcast(data=stool.result, Cytokine ~ Richness, value.var =  "beta")
row.names(wide.stool) <- wide.stool$Cytokine
wide.stool <- select(wide.stool, -Cytokine)

wide.ptable.stool <- dcast(data=stool.result, Cytokine ~ Richness, value.var =  "adj.P")
row.names(wide.ptable.stool) <- wide.ptable.stool$Cytokine
wide.ptable.stool <- select(wide.ptable.stool, -Cytokine)

#draw heatmap
col_fun = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))

row_dend = as.dendrogram(hclust(dist(wide.stool)))
row_dend = color_branches(row_dend, k = 4)

heatmap.stool <- Heatmap(wide.stool, col = col_fun,
                         column_title = "Stool_Richness:Cytokine",
                         column_km = 2,
                         column_dend_side = "top",
                         cluster_rows = row_dend, 
                         row_dend_reorder = TRUE,
                         row_names_gp=grid::gpar(fontsize = 9,fontface = "bold"), 
                         column_names_gp = grid::gpar(fontsize = 10),
                         cell_fun = function(j, i, x, y, w, h, fill){
                           gb = textGrob("*")
                           gb_w = convertWidth(grobWidth(gb), "mm")
                           gb_h = convertHeight(grobHeight(gb), "mm")
                           if(wide.ptable.stool[i, j] < 0.2) {
                             grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)}
                           })
heatmap.stool

pdf("./Stool.richness.cytokine.pdf", height = 10, width = 4.5)
draw(heatmap.stool)
dev.off()

#skin
wide.skin <- dcast(data=skin.result, Cytokine ~ Richness, value.var =  "beta")
row.names(wide.skin) <- wide.skin$Cytokine
wide.skin <- select(wide.skin, -Cytokine)

wide.ptable.skin <- dcast(data=skin.result, Cytokine ~ Richness, value.var =  "adj.P")
row.names(wide.ptable.skin) <- wide.ptable.skin$Cytokine
wide.ptable.skin <- select(wide.ptable.skin, -Cytokine)

#draw heatmap
col_fun = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))

row_dend = as.dendrogram(hclust(dist(wide.skin)))
row_dend = color_branches(row_dend, k = 5)

heatmap.skin <- Heatmap(wide.skin, col = col_fun,
                         column_title = "Skin_Richness:Cytokine",
                         column_km = 2,
                         column_dend_side = "top",
                         cluster_rows = row_dend, 
                         row_dend_reorder = TRUE,
                         row_names_gp=grid::gpar(fontsize = 9,fontface = "bold"), 
                         column_names_gp = grid::gpar(fontsize = 10),
                         cell_fun = function(j, i, x, y, w, h, fill){
                           gb = textGrob("*")
                           gb_w = convertWidth(grobWidth(gb), "mm")
                           gb_h = convertHeight(grobHeight(gb), "mm")
                           if(wide.ptable.skin[i, j] < 0.2) {
                             grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)}
                         })
heatmap.skin

pdf("./skin.richness.cytokine.pdf", height = 10, width = 4.5)
draw(heatmap.skin)
dev.off()

#oral
wide.oral <- dcast(data=oral.result, Cytokine ~ Richness, value.var =  "beta")
row.names(wide.oral) <- wide.oral$Cytokine
wide.oral <- select(wide.oral, -Cytokine)

wide.ptable.oral <- dcast(data=oral.result, Cytokine ~ Richness, value.var =  "adj.P")
row.names(wide.ptable.oral) <- wide.ptable.oral$Cytokine
wide.ptable.oral <- select(wide.ptable.oral, -Cytokine)

#draw heatmap
col_fun = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))

row_dend = as.dendrogram(hclust(dist(wide.oral)))
row_dend = color_branches(row_dend, k = 4)

heatmap.oral <- Heatmap(wide.oral, col = col_fun,
                        column_title = "Oral_Richness:Cytokine",
                        column_km = 1,
                        column_dend_side = "top",
                        cluster_rows = row_dend, 
                        row_dend_reorder = TRUE,
                        row_names_gp=grid::gpar(fontsize = 9,fontface = "bold"), 
                        column_names_gp = grid::gpar(fontsize = 10),
                        cell_fun = function(j, i, x, y, w, h, fill){
                          gb = textGrob("*")
                          gb_w = convertWidth(grobWidth(gb), "mm")
                          gb_h = convertHeight(grobHeight(gb), "mm")
                          if(wide.ptable.oral[i, j] < 0.2) {
                            grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)}
                        })
heatmap.oral

pdf("./oral.richness.cytokine.pdf", height = 10, width = 4)
draw(heatmap.oral)
dev.off()

#nasal
wide.nasal <- dcast(data=nasal.result, Cytokine ~ Richness, value.var =  "beta")
row.names(wide.nasal) <- wide.nasal$Cytokine
wide.nasal <- select(wide.nasal, -Cytokine)

wide.ptable.nasal <- dcast(data=nasal.result, Cytokine ~ Richness, value.var =  "adj.P")
row.names(wide.ptable.nasal) <- wide.ptable.nasal$Cytokine
wide.ptable.nasal <- select(wide.ptable.nasal, -Cytokine)

#draw heatmap
col_fun = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))

row_dend = as.dendrogram(hclust(dist(wide.nasal)))
row_dend = color_branches(row_dend, k = 4)

heatmap.nasal <- Heatmap(wide.nasal, col = col_fun,
                        column_title = "nasal_Richness:Cytokine",
                        column_km = 2,
                        column_dend_side = "top",
                        cluster_rows = row_dend, 
                        row_dend_reorder = TRUE,
                        row_names_gp=grid::gpar(fontsize = 9,fontface = "bold"), 
                        column_names_gp = grid::gpar(fontsize = 10),
                        cell_fun = function(j, i, x, y, w, h, fill){
                          gb = textGrob("*")
                          gb_w = convertWidth(grobWidth(gb), "mm")
                          gb_h = convertHeight(grobHeight(gb), "mm")
                          if(wide.ptable.nasal[i, j] < 0.2) {
                            grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)}
                        })
heatmap.nasal

pdf("./nasal.richness.cytokine.pdf", height = 10, width = 4)
draw(heatmap.nasal)
dev.off()





