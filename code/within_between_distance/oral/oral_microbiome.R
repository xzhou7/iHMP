###
no_function()

setwd(masstools::get_project_wd())
library(tidyverse)
library(phyloseq)
rm(list = ls())
####load raw data
load("data/from_xin/DetailedPhyloseq.RData")
ls()

dir.create("data_analysis/oral_microbiome/")
dir.create("data_analysis/oral_microbiome/cladogram/")
setwd("data_analysis/oral_microbiome/cladogram/")

load("bc_dist_by_genus_OR.rda")
bc_dist_oral = bc_dist_by_genus
rm(list = "bc_dist_by_genus")
dim(bc_dist_oral)
unique(bc_dist_oral$genus)

###remove distance == 1
bc_dist_oral =
  bc_dist_oral %>%
  dplyr::filter(bc_dist != 1)

dim(bc_dist_oral)

##remove genus with between or within distance sample number < 3
library(plyr)

number =
  bc_dist_oral %>%
  dplyr::group_by(genus, type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  plyr::dlply(.variables = .(genus)) %>%
  purrr::map(function(x) {
    x =
      x %>%
      dplyr::filter(n >= 10)
    
    if (all(x$type != "within") | all(x$type != "between")) {
      return(NULL)
    }
    
    paste(x$genus, x$type, sep = "_")
  }) %>%
  unlist() %>%
  unique()

bc_dist_oral2 = 
  bc_dist_oral %>% 
  dplyr::mutate(name = paste(genus, type , sep = "_")) %>% 
  dplyr::filter(name %in% number)

dim(bc_dist_oral2)

##calculate the median distance for each genuse
temp_data =
  bc_dist_oral2 %>%
  dplyr::group_by(genus, type) %>%
  dplyr::summarise(median = median(bc_dist),
                   sd = sd(bc_dist)) 

target_name = 
  temp_data %>% 
  dplyr::arrange(desc(median)) %>% 
  dplyr::filter(type == "within") %>% 
  dplyr::pull(genus) %>% 
  head(20)

#####boxplot
bc_dist_oral2 %>% 
  dplyr::filter(genus %in% target_name) %>%
  dplyr::mutate(genus = factor(genus, levels = rev(target_name))) %>% 
  ggplot(aes(x = bc_dist, y = genus)) +
  geom_boxplot(aes(color = type), 
               outlier.shape = NA)

# temp_data %>% 
#   dplyr::filter(genus == "Peptostreptococcus")

library(plyr)

####here we change fc1 as between - minus
temp_data = 
  temp_data %>%
  plyr::dlply(.variables = "genus") %>%
  purrr::map(function(x) {
    cat(x$genus[1], " ")
    if(nrow(x) == 1){
      return(NULL)
    }
    between = x$median[x$type == "between"] 
    within = x$median[x$type == "within"]
    family = x$median[x$type == "family"]
    fc1 = between - within
    fc2 = (family - within) / (between - family)
    if(length(fc2) == 0){
      fc2 = NA
    }else{
      if(fc2 < 0){
        fc2 = 0
      }
      if(fc2 > 1){
        fc2 = 1
      }  
    }
    c(fc1 = fc1, fc2 = fc2)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

# temp_data %>% 
#   dplyr::arrange(desc(fc1)) %>% 
#   head()

####get fc1_p
fc1_p =
  bc_dist_oral2 %>%
  dplyr::filter(type != "family") %>%
  plyr::dlply(.variables = .(genus)) %>%
  purrr::map(function(x) {
    cat(x$genus[1], " ")
    p = 
      tryCatch(wilcox.test(x$bc_dist[x$type == "between"],
                           x$bc_dist[x$type == "within"])$p.value, error = function(x){
                             NA
                           })
    p
  }) %>% 
  unlist() %>% 
  data.frame(fc1_p = .) %>% 
  tibble::rownames_to_column(var = "Genus") %>% 
  dplyr::filter(!is.na(fc1_p))

fc1_p$fc1_p_adjust = p.adjust(fc1_p$fc1_p, method = "BH")

temp_data =
  temp_data %>% 
  rownames_to_column(var = "Genus")  %>% 
  dplyr::left_join(fc1_p, by = "Genus")

####HMP
physeqGenus_OR

expression_data = physeqGenus_OR@otu_table@.Data %>%
  t() %>%
  as.data.frame()

variable_info = as.data.frame(physeqGenus_OR@tax_table@.Data)

sample_info = get_variable(physeq = physeqGenus_OR)

####remove some variables
###change to relative counts
##only remain taxs with median > 1e-5
# physeqGenus_OR2 = filter_taxa(physeq = physeqGenus_OR, flist = function(x) median(x) > 1e-5, prune =  TRUE)
##remove some variables which genus are not in temp_data
# physeqGenus_OR2 = filter_taxa(physeq = physeqGenus_OR, 
#                               flist = function(x) median(x) > 1e-5, prune =  TRUE)
#                               
# physeqGenus_OR2 =
#   subset_taxa(physeq = physeqGenus_OR, Genus %in% grep("Unclassified",temp_data$Genus, value = TRUE, invert = TRUE))

physeqGenus_OR2 =
  subset_taxa(physeq = physeqGenus_OR, Genus %in% temp_data$Genus)

# ##remove some genus which have no root
# ##
# physeqGenus_OR2 =
#   subset_taxa(physeq = physeqGenus_OR2, !Genus %in% "Unclassified_Actinobacteria")

library(microbiomeViz)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)

GP <- physeqGenus_OR2
GP
GP = fix_duplicate_tax(GP)
GP

###create tree
tree = microbiomeViz::parsePhyloseq(
  physeq = GP,
  use_abundance = TRUE,
  node.size.scale = 1,
  node.size.offset = 1
)

raw_p <- tree.backbone(
  tree = tree,
  size = 0.5,
  shape = 21,
  layout = "circular",
  fill = "white",
  color = "black"
)

raw_p

###add tip label (genus label)
raw_p$data$label2 = raw_p$data$label %>% 
  stringr::str_split("__") %>% 
  purrr::map(function(x){
    x[2]
  }) %>% 
  unlist()

raw_p$data$label2[as.character(raw_p$data$nodeClass) != "g"] = NA

new_info =
  data.frame(Genus = raw_p$data$label2) %>%
  dplyr::left_join(variable_info[, c("Genus", "Kingdom", "Phylum", "Class", "Order", "Family")],
                   by = "Genus") %>%
  dplyr::left_join(temp_data, by = "Genus") 

raw_p$data =
  cbind(raw_p$data, new_info)

raw_p$data =
  raw_p$data %>%
  dplyr::mutate(
    fc1_star =
      case_when(
        is.na(fc1_p_adjust) ~ "",
        fc1_p_adjust > 0.05 ~ "",
        fc1_p_adjust <= 0.05 &
          fc1_p_adjust > 0.01 ~ "*",
        fc1_p_adjust <= 0.01 &
          fc1_p_adjust > 0.001 ~ "*",
        fc1_p_adjust <= 0.001 ~ "*"
      )
  )

# p$data$fc1[which(p$data$fc1_p_adjust > 0.05)] = NA

##add tip lab
##only add some tip points
idx1 = which(raw_p$data$fc1_p_adjust > 0.05)

raw_p$data$label2[idx1] = NA

p1 =
  raw_p +
  geom_tiplab(
    aes(label = label2,
        color = Phylum),
    offset = 1.3,
    size = 2,
    show.legend = FALSE
  )

p1

######add heatmap
##add fc1 information
p1$data$fc1[which(p1$data$fc1 < 0)] = 0

##add multiple tip information
fc1_info =
  p1$data[, c("fc1", "isTip")]

rownames(fc1_info) = p1$data$label

fc1_info = fc1_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

fc1_info$fc1[which(fc1_info$fc1 < 0)] = 0

p2 =
  gheatmap(
    p = p1,
    data = fc1_info,
    offset = -0.1,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25, 
    colnames = FALSE,
    color = "black"
  ) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tippoint(
    mapping = aes(shape = fc1_star),
    x = 6.35,
    size = 3,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = c("*" = "*")) 

p2

##add fc2 information
##add multiple tip information
fc2_info =
  p2$data[, c("fc2", "isTip")] 
  # dplyr::rename(fc2 = fc1)

rownames(fc2_info) = p2$data$label

fc2_info = fc2_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

fc2_info$fc2[which(fc2_info$fc2 < 0)] = 0

# idx = sample(1:150, 100)
# fc2_info$fc2[idx] = -fc2_info$fc2[idx]

p3 = p2 +
  ggnewscale::new_scale_fill()

p3 =
  gheatmap(
    p = p3,
    data = fc2_info,
    offset = 0.5,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25, 
    colnames = FALSE, color = "grey"
  ) +
  scale_fill_gradient2(low = ggsci::pal_aaas()(n=10)[3], 
                       mid = "white", 
                       high = ggsci::pal_aaas()(n=10)[4]) 

p3

plot(p3$data$fc1_p_adjust)

####output results
temp_data = 
  p3$data

temp_data = 
  temp_data %>%
  dplyr::filter(!is.na(fc1_p_adjust))

openxlsx::write.xlsx(temp_data,
                     file = "tree_data.xlsx",
                     asTable = TRUE,
                     overwrite = TRUE)

table(temp_data$Kingdom)
table(temp_data$Phylum)
table(temp_data$Class)
table(temp_data$Order)
table(temp_data$Family)

###Kingdom level
Kingdom = unique(temp_data$Kingdom)
if(length(Kingdom) == 1){
  Kingdom_result =
    data.frame(
      class1 = NA,
      class2 = NA,
      class1_sig_ratio = NA,
      class2_sig_ratio = NA,
      p = NA
    )
}else{
  Kingdom_result =
    purrr::map(1:(length(Kingdom)-1), function(i) {
      cat(i, "")
      purrr::map((i+1):length(Kingdom), .f = function(j){
        data1 =
          temp_data[temp_data$Kingdom == Kingdom[i],]
        
        data2 =
          temp_data[temp_data$Kingdom == Kingdom[j],]
        
        test = rbind(
          c(sum(data1$fc1_p_adjust < 0.05),
            sum(data1$fc1_p_adjust >= 0.05)),
          c(sum(data2$fc1_p_adjust < 0.05),
            sum(data2$fc1_p_adjust >= 0.05))
        )
        
        class1_sig_ratio = test[1,1]/sum(test[1,])
        class2_sig_ratio = test[2,1]/sum(test[2,])
        
        if(nrow(data1) < 10 | nrow(data2) < 10){
          return(data.frame(class1 = Kingdom[i],
                            class2 = Kingdom[j],
                            class1_sig_ratio = NA,
                            class2_sig_ratio = NA,
                            p = NA))
        }
        
        p =
          chisq.test(test)$p.value
        
        data.frame(class1 = Kingdom[i],
                   class2 = Kingdom[j],
                   class1_sig_ratio = class1_sig_ratio,
                   class2_sig_ratio = class1_sig_ratio,
                   p = NA)
      }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    }) %>%
    as.data.frame()  
}

Kingdom_result

###Phylum level
Phylum = unique(temp_data$Phylum)
Phylum_result =
purrr::map(1:(length(Phylum)-1), function(i) {
  cat(i, "")
  purrr::map((i+1):length(Phylum), .f = function(j){
    data1 =
      temp_data[temp_data$Phylum == Phylum[i],]

    data2 =
      temp_data[temp_data$Phylum == Phylum[j],]

    test = rbind(
      c(sum(data1$fc1_p_adjust < 0.05),
        sum(data1$fc1_p_adjust >= 0.05)),
      c(sum(data2$fc1_p_adjust < 0.05),
        sum(data2$fc1_p_adjust >= 0.05))
    )

    class1_sig_ratio = test[1,1]/sum(test[1,])
    class2_sig_ratio = test[2,1]/sum(test[2,])

    if(nrow(data1) < 10 | nrow(data2) < 10){
      return(data.frame(class1 = Phylum[i],
                        class2 = Phylum[j],
                        class1_sig_ratio = NA,
                        class2_sig_ratio = NA,
                        p = NA))
    }

    p =
      chisq.test(test)$p.value

    data.frame(class1 = Phylum[i],
               class2 = Phylum[j],
               class1_sig_ratio = class1_sig_ratio,
               class2_sig_ratio = class2_sig_ratio,
               p = p)
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
}) %>%
  do.call(rbind, .) %>%
  as.data.frame()

Phylum_result %>%
  dplyr::filter(!is.na(p))

###Class level
Class = unique(temp_data$Class)
table(temp_data$Class)
Class_result =
  purrr::map(1:(length(Class)-1), function(i) {
    cat(i, "")
    purrr::map((i+1):length(Class), .f = function(j){
      data1 =
        temp_data[temp_data$Class == Class[i],]

      data2 =
        temp_data[temp_data$Class == Class[j],]

      test = rbind(
        c(sum(data1$fc1_p_adjust < 0.05),
          sum(data1$fc1_p_adjust >= 0.05)),
        c(sum(data2$fc1_p_adjust < 0.05),
          sum(data2$fc1_p_adjust >= 0.05))
      )

      class1_sig_ratio = test[1,1]/sum(test[1,])
      class2_sig_ratio = test[2,1]/sum(test[2,])

      if(nrow(data1) < 10 | nrow(data2) < 10){
        return(data.frame(class1 = Class[i],
                          class2 = Class[j],
                          class1_sig_ratio = NA,
                          class2_sig_ratio = NA,
                          p = NA))
      }

      p =
        chisq.test(test)$p.value

      data.frame(class1 = Class[i],
                 class2 = Class[j],
                 class1_sig_ratio = class1_sig_ratio,
                 class2_sig_ratio = class2_sig_ratio,
                 p = p)
    }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

Class_result %>%
  dplyr::filter(!is.na(p))

###Order level
Order = unique(temp_data$Order)
table(temp_data$Order)
Order_result =
  purrr::map(1:(length(Order)-1), function(i) {
    cat(i, "")
    purrr::map((i+1):length(Order), .f = function(j){
      data1 =
        temp_data[temp_data$Order == Order[i],]

      data2 =
        temp_data[temp_data$Order == Order[j],]

      test = rbind(
        c(sum(data1$fc1_p_adjust < 0.05),
          sum(data1$fc1_p_adjust >= 0.05)),
        c(sum(data2$fc1_p_adjust < 0.05),
          sum(data2$fc1_p_adjust >= 0.05))
      )

      class1_sig_ratio = test[1,1]/sum(test[1,])
      class2_sig_ratio = test[2,1]/sum(test[2,])

      if(nrow(data1) < 10 | nrow(data2) < 10){
        return(data.frame(class1 = Order[i],
                          class2 = Order[j],
                          class1_sig_ratio = NA,
                          class2_sig_ratio = NA,
                          p = NA))
      }

      p =
        chisq.test(test)$p.value

      data.frame(
        class1 = Order[i],
        class2 = Order[j],
        class1_sig_ratio = class1_sig_ratio,
        class2_sig_ratio = class2_sig_ratio,
        p = p
      )
    }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

Order_result %>%
  dplyr::filter(!is.na(p))

###Family level
Family = unique(temp_data$Family)
table(temp_data$Family)
Family_result =
  purrr::map(1:(length(Family)-1), function(i) {
    cat(i, "")
    purrr::map((i+1):length(Family), .f = function(j){
      data1 =
        temp_data[temp_data$Family == Family[i],]

      data2 =
        temp_data[temp_data$Family == Family[j],]


      test = rbind(
        c(sum(data1$fc1_p_adjust < 0.05),
          sum(data1$fc1_p_adjust >= 0.05)),
        c(sum(data2$fc1_p_adjust < 0.05),
          sum(data2$fc1_p_adjust >= 0.05))
      )

      class1_sig_ratio = test[1,1]/sum(test[1,])
      class2_sig_ratio = test[2,1]/sum(test[2,])

      if(nrow(data1) < 10 | nrow(data2) < 10){
        return(data.frame(class1 = Family[i],
                          class2 = Family[j],
                          class1_sig_ratio = NA,
                          class2_sig_ratio = NA,
                          p = NA))
      }

      p =
        chisq.test(test)$p.value

      data.frame(
        class1 = Family[i],
        class2 = Family[j],
        class1_sig_ratio = class1_sig_ratio,
        class2_sig_ratio = class2_sig_ratio,
        p = p
      )

    }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

######class
save(Kingdom_result, file = "Kingdom_result")
save(Phylum_result, file = "Phylum_result")
save(Class_result, file = "Class_result")
save(Order_result, file = "Order_result")
save(Family_result, file = "Family_result")

Kingdom_result =
Kingdom_result %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::mutate(p_adjust = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p_adjust)

Phylum_result =
Phylum_result %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::mutate(p_adjust = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p_adjust)

Class_result =
Class_result %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::mutate(p_adjust = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p_adjust)

Order_result =
Order_result %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::mutate(p_adjust = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p_adjust)

Family_result =
Family_result %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::mutate(p_adjust = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p_adjust)

####output result
library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")

addWorksheet(wb, sheetName = "Kingdom", gridLines = TRUE)
addWorksheet(wb, sheetName = "Phylum", gridLines = TRUE)
addWorksheet(wb, sheetName = "Class", gridLines = TRUE)
addWorksheet(wb, sheetName = "Order", gridLines = TRUE)
addWorksheet(wb, sheetName = "Family", gridLines = TRUE)

freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
freezePane(wb, sheet = 3, firstRow = TRUE, firstCol = TRUE)
freezePane(wb, sheet = 4, firstRow = TRUE, firstCol = TRUE)
freezePane(wb, sheet = 5, firstRow = TRUE, firstCol = TRUE)

writeDataTable(wb, sheet = 1, x = Kingdom_result, colNames = TRUE,rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = Phylum_result, colNames = TRUE,rowNames = FALSE)
writeDataTable(wb, sheet = 3, x = Class_result, colNames = TRUE,rowNames = FALSE)
writeDataTable(wb, sheet = 4, x = Order_result, colNames = TRUE,rowNames = FALSE)
writeDataTable(wb, sheet = 5, x = Family_result, colNames = TRUE,rowNames = FALSE)

saveWorkbook(wb, "chisq.test.result.xlsx", overwrite = TRUE)

Kingdom_result
Phylum_result
Class_result
Order_result
Family_result

###highlight
final_p =
  p3 
  # geom_hilight(
  #   node = p3$data$node[which(p3$data$label == "c__Bacteroidia")],
  #   fill = ggsci::pal_lancet()(n = 9)[1],
  #   alpha = .4
  # ) +
  # geom_hilight(
  #   node = p3$data$node[which(p3$data$label == "c__Clostridia")],
  #   fill = ggsci::pal_lancet()(n = 9)[3],
  #   alpha = .4
  # ) +
  # geom_hilight(
  #   node = p3$data$node[which(p3$data$label == "c__Betaproteobacteria")],
  #   fill = ggsci::pal_lancet()(n = 9)[5],
  #   alpha = .4
  # )

final_p

ggsave(final_p, filename = "oral.tree.pdf", width = 12, height = 12)
