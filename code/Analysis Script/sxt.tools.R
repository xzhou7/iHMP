iris_color = 
  c("IR" = ggsci::pal_nejm()(n=10)[2],
    "IS" = ggsci::pal_nejm()(n=10)[1]) 


season_color = c(
  "spring" = ggsci::pal_locuszoom()(n=7)[3],
  "summer" = ggsci::pal_locuszoom()(n=7)[1],
  "autumn" = ggsci::pal_locuszoom()(n=7)[2],
  "winter" = ggsci::pal_locuszoom()(n=7)[4]
)


omics_color = c(
  "Metabolite" = ggsci::pal_jama()(n=7)[1],
  "Lipidome" = ggsci::pal_aaas()(n=7)[2],
  "Stool microbiome" = ggsci::pal_jama()(n=7)[2],
  "Skin microbiome" = ggsci::pal_jama()(n=7)[3],
  "Oral microbiome" = ggsci::pal_jama()(n=7)[4],
  "Nasal microbiome" = ggsci::pal_jama()(n=7)[5],
  "Cytokine" = ggsci::pal_jama()(n=7)[6],
  "Exposome" = ggsci::pal_jama()(n=7)[7],
  "Proteome" = ggsci::pal_d3()(n=7)[2]
)

omics_shape = c(
  "Metabolite" = 6,
  "Lipidome" = 11,
  "Stool microbiome" = 15,
  "Skin microbiome" = 15,
  "Oral microbiome" = 15,
  "Nasal microbiome" = 15,
  "Cytokine" = 5,
  "Exposome" = 28,
  "Proteome" = 13
)


metabolite_class_color = c(
  "Amino Acid" = ggsci::pal_lancet()(n=8)[1],
  "Carbohydrate" = ggsci::pal_lancet()(n=8)[2],
  "Cofactors and Vitamins" = ggsci::pal_lancet()(n=8)[3],
  "Energy" = ggsci::pal_lancet()(n=8)[4],
  "Lipid" = ggsci::pal_lancet()(n=8)[5],
  "Nucleotide" = ggsci::pal_lancet()(n=8)[6],
  "Peptide" = ggsci::pal_lancet()(n=8)[7],
  "Xenobiotics" = ggsci::pal_lancet()(n=8)[8]
)

lipid_class_color =
  c(
    "CE" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[1],
    "CER" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[2],
    "DAG" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[3],
    "FFA" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[4],
    "LPC" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[5],
    "LPE" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[6],
    "PC" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[7],
    "PE" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[8],
    "PI" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[9],
    "SM" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[10],
    "TAG" = RColorBrewer::brewer.pal(n = 11, name = "Set3")[11]
  )


between_within_color = 
  c("between" = ggsci::pal_nejm()(n=8)[1],
    "within" = ggsci::pal_nejm()(n=8)[2],
    "family" = ggsci::pal_nejm()(n=8)[3])


base_theme =
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())


body_site_color = c(
  "Stool" = ggsci::pal_jama()(n=7)[2],
  "Skin" = ggsci::pal_jama()(n=7)[3],
  "Oral" = ggsci::pal_jama()(n=7)[4],
  "Nasal" = ggsci::pal_jama()(n=7)[5]
)


infection_color = 
  c("Healthy" = ggsci::pal_npg()(n=7)[3],
    "Infection" = ggsci::pal_npg()(n=7)[1])

phylum_name =
  c(
    "Actinobacteria",
    "Bacteroidetes",
    "Cyanobacteria/Chloroplast"  ,
    "Firmicutes",
    "Lentisphaerae",
    "Proteobacteria",
    "Synergistetes",
    "Verrucomicrobia",
    "Campilobacterota",
    "Candidatus_Saccharibacteria",
    "Fusobacteria",
    "Plantae",
    "Tenericutes",
    "Spirochaetes"
  )

phylum_color = 
  ggsci::pal_simpsons()(n=length(phylum_name))

names(phylum_color) = phylum_name

plot_pvca <- function(object){
  data <- data.frame(label = object$label, 
                     value = object$dat[1,] * 100,
                     stringsAsFactors = FALSE) %>% 
    dplyr::arrange(value) %>% 
    dplyr::mutate(label = factor(label, levels = label)) %>% 
    dplyr::mutate(
      class = case_when(
        stringr::str_detect(label, ":") ~ "no",
        TRUE ~ "yes"
      )
    )
  
  plot <- 
    data %>% 
    ggplot(aes(value, label)) +
    geom_segment(aes(x = 0, y = label, xend = value, yend = label)) +
    geom_point(size = 4,
               shape = 21, 
               aes(fill = class),
               show.legend = FALSE) +
    scale_fill_manual(values = c(
      "yes" = ggsci::pal_aaas()(n=10)[2],
      "no" = ggsci::pal_aaas()(n=10)[1]
    )) +
    geom_text(aes(x = value, y = label, 
                  label = round(value,2)), 
              hjust = -0.3, size = 4) +
    labs(x = "Weighted average proportion variance (%)",
         y = "") +
    scale_x_continuous(expand = expansion(mult = c(0.05,0.2))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12))
  
  plot  
}



#' A function for principal variance component analysis
#'
#' The function is written based on the 'pvcaBatchAssess' function of the PVCA R package
#' and slightly changed to make it more efficient and flexible for sequencing read counts data.
#' (http://watson.nci.nih.gov/bioc_mirror/packages/release/bioc/manuals/pvca/man/pvca.pdf)
#'
#' @param counts The Normalized(e.g. TMM)/ log-transformed reads count matrix from sequencing data (row:gene/feature, col:sample)
#' @param meta  The Meta data matrix containing predictor variables (row:sample, col:predictor)
#' @param threshold The proportion of the variation in read counts explained by top k PCs. This value determines the number of top PCs to be used in pvca.
#' @param inter TRUE/FALSE - include/do not include pairwise interactions of predictors
#'
#' @return std.prop.val The vector of proportions of variation explained by each predictor.
#'
#' @export
#'

PVCA <- function(counts, meta, threshold, inter){
  library(nlme)
  library(lme4)
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)
  
  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}
  
  pred.list <- colnames(meta)
  meta <- droplevels(meta)
  
  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
  
  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit)
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}


PlotPVCA <- function(pvca.res, title){
  plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res * 100)
  plot.dat <- 
    plot.dat %>% 
    dplyr::arrange(prop) %>% 
    dplyr::mutate(eff = factor(eff, levels = eff))
  
  plot.dat %>% 
    ggplot(aes(prop, eff)) +
    geom_bar(stat = "identity", fill = ggsci::pal_aaas()(n=10)[3],
             color = "black") +
    geom_text(aes(x = prop + 1.5, y = eff, label = round(prop,2))) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_discrete(breaks = c("mother_pre_bmi", "resid", "alcohol_during_preg_yes_no", 
                                "child_weight_g", "mother_age", "delivery_ga", 
                                "parity", "child_length_cm", "g_stage",
                                "child_sex", "season", "days", "subject_id"),
                     labels = c("Mother BMI", "Residual", "Alcohol",
                                "Child weight", "Mother age", "Delivery GA",
                                "Parity", "Child length", "GA",
                                "Child sex", "Begin season", "Days", "Subject ID")) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 13)
    ) +
    labs(x = "Weighted average proportion variance (%)",
         y = "")
}


####calculate the liner mixed model adjusted spearman correlation for microbiome and metabolome
lm_adjusted_cor = function(data_set1,
                           data_set2,
                           sample_info,
                           method = c("spearman", "pearson", "all"),
                           threads = 5) {
  method = match.arg(method)
  library(future)
  library(furrr)
  # plan(strategy = multisession(workers = threads))
  
  cor =
    rownames(data_set1) %>%
    purrr::map(function(name) {
      cat(name, " ")
      x = as.numeric(data_set1[name, ])
      
      temp_cor =
        purrr::map(
          as.data.frame(t(data_set2)),
          .f = function(y) {
            temp_data =
              data.frame(x = x,
                         y = y,
                         sample_info)
            temp_data$Gender[temp_data$Gender == 'F'] = 0
            temp_data$Gender[temp_data$Gender == 'M'] = 1
            temp_data$Gender = as.numeric(temp_data$Gender)
            
            temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
            temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
            temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
            temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
            temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
            
            temp_data$SSPG = as.numeric(temp_data$SSPG)
            temp_data$FPG = as.numeric(temp_data$FPG)
            
            ##only remain the subjects with less than 5 time points
            remain_subject_id =
              temp_data %>%
              dplyr::select(subject_id) %>%
              dplyr::group_by(subject_id) %>%
              dplyr::summarise(n = n()) %>%
              dplyr::filter(n >= 5) %>%
              pull(subject_id)
            
            temp_data =
              temp_data %>%
              dplyr::filter(subject_id %in% remain_subject_id)
            
            ##linear mixed model to adjusted the x and y
            adjusted_x =
              lme4::lmer(formula = x ~ Gender + Adj.age + Ethnicity + (1 |
                                                                         subject_id),
                         data = temp_data) %>%
              residuals()
            
            adjusted_y =
              lme4::lmer(
                formula = y ~ Gender + Adj.age + Ethnicity + (1 |
                                                                subject_id),
                data = temp_data
              ) %>%
              residuals()
            
            if (method == "all") {
              cor_value1 =
                cor.test(adjusted_x, adjusted_y, method = "spearman")
              
              cor_value2 =
                cor.test(adjusted_x, adjusted_y, method = "pearson")
              
              result1 =
                c(cor = unname(cor_value1$estimate),
                  p = unname(cor_value1$p.value))
              
              if (is.na(result1[1])) {
                result1[1] = 0
              }
              
              if (is.na(result1[2])) {
                result1[2] = 1
              }
              
              result2 =
                c(cor = unname(cor_value2$estimate),
                  p = unname(cor_value2$p.value))
              
              if (is.na(result2[1])) {
                result2[1] = 0
              }
              
              if (is.na(result2[2])) {
                result2[2] = 1
              }
              
              list(result1 = result1, result2 = result2)
              
            } else {
              cor_value =
                cor.test(adjusted_x, adjusted_y, method = method)
              
              result =
                c(cor = unname(cor_value$estimate),
                  p = unname(cor_value$p.value))
              
              if (is.na(result[1])) {
                result[1] = 0
              }
              
              if (is.na(result[2])) {
                result[2] = 1
              }
              list(result)
            }
          }
        ) 
      # do.call(rbind, .) %>%
      # as.data.frame()
      
      if(method == "all") {
        temp_cor1 = 
          temp_cor %>% 
          purrr::map(function(x){
            x[[1]]
          }) %>% 
          do.call(rbind, .) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column(var = "metabolite") %>%
          dplyr::mutate(microbiome = name) %>%
          dplyr::select(microbiome, metabolite, dplyr::everything())
        
        temp_cor1$p_adjust = p.adjust(temp_cor1$p, method = "BH")
        
        temp_cor2 = 
          temp_cor %>% 
          purrr::map(function(x){
            x[[2]]
          }) %>% 
          do.call(rbind, .) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column(var = "metabolite") %>%
          dplyr::mutate(microbiome = name) %>%
          dplyr::select(microbiome, metabolite, dplyr::everything())
        
        temp_cor2$p_adjust = p.adjust(temp_cor2$p, method = "BH")
        
        list(temp_cor1 = temp_cor1,
             temp_cor2 = temp_cor2)
      }else{
        temp_cor =
          temp_cor %>%
          tibble::rownames_to_column(var = "metabolite") %>%
          dplyr::mutate(microbiome = name) %>%
          dplyr::select(microbiome, metabolite, dplyr::everything())
        temp_cor$p_adjust = p.adjust(temp_cor$p, method = "BH")
        list(temp_cor)
      }
    })
  # do.call(rbind, .) %>%
  # as.data.frame()
  
  if(method == "all"){
    cor1 = 
      cor %>% 
      purrr::map(function(x){
        x[[1]]
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    cor2 = 
      cor %>% 
      purrr::map(function(x){
        x[[2]]
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    cor1$p_adjust2 =
      p.adjust(cor1$p, method = "BH")
    
    cor2$p_adjust2 =
      p.adjust(cor2$p, method = "BH")
    list(spearman = cor1, 
         pearson = cor2)
  }else{
    cor$p_adjust2 =
      p.adjust(cor$p, method = "BH")
    list(cor)
  }
}


####calculate the liner mixed model adjusted spearman correlation intra omics data
intra_lm_adjusted_cor = 
  function(expression_data,
           sample_info) {
    
  library(future)
  library(furrr)
  plan(strategy = multisession(workers = 5))
  
  cor = 
    1:(nrow(expression_data)-1) %>% 
    furrr::future_map(function(i){
      # cat(i, " ")
      x = as.numeric(expression_data[i, ])
      
      temp_cor = 
      purrr::map(
        i:nrow(expression_data),
        .f = function(j) {
          y = as.numeric(expression_data[j, ])
          temp_data =
            data.frame(x = x,
                       y = y,
                       sample_info)
          
          temp_data$Gender[temp_data$Gender == 'F'] = 0
          temp_data$Gender[temp_data$Gender == 'M'] = 1
          temp_data$Gender = as.numeric(temp_data$Gender)
          
          temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
          temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
          temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
          temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
          temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
          
          temp_data$SSPG = as.numeric(temp_data$SSPG)
          temp_data$FPG = as.numeric(temp_data$FPG)
          
          ##only remain the subjects with less than 5 time points
          remain_subject_id =
            temp_data %>%
            dplyr::select(subject_id) %>%
            dplyr::group_by(subject_id) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::filter(n >= 5) %>%
            pull(subject_id)
          
          temp_data =
            temp_data %>%
            dplyr::filter(subject_id %in% remain_subject_id)
          
          ##linear mixed model to adjusted the x and y
          adjusted_x =
            lme4::lmer(
              formula = x ~ Gender + Adj.age + Ethnicity + (1 |
                                                              subject_id),
              data = temp_data
            ) %>%
            residuals()
          
          adjusted_y =
            lme4::lmer(
              formula = y ~ Gender + Adj.age + Ethnicity + (1 |
                                                              subject_id),
              data = temp_data
            ) %>%
            residuals()
          
          cor_value =
            cor.test(adjusted_x, adjusted_y, method = "spearman")
          
          result =
            data.frame(
              name1 = rownames(expression_data)[i],
              name2 = rownames(expression_data)[2],
              cor = unname(cor_value$estimate),
              p = unname(cor_value$p.value)
            )
          
          if (is.na(result$cor)) {
            result$cor = 0
          }
          
          if (is.na(result$p)) {
            result$p = 1
          }
          result
        }
      ) %>% 
        do.call(rbind, .) %>% 
        as.data.frame()
      
        
      temp_cor$p_adjust = p.adjust(temp_cor$p, method = "BH")
      temp_cor
    }, .progress = TRUE) %>% 
    do.call(rbind, .) %>%
    as.data.frame()
  
  cor$p_adjust2 =
    p.adjust(cor$p, method = "BH")
  cor
}


####calculate the partial spearman correlation for microbiome and metabolome
partial_cor = function(microbiome_data,
                       metabolome_data,
                       sample_info,
                       method = c("spearman", "pearson", "all"),
                       threads = 5) {
  # browser()
  method = match.arg(method)
  library(future)
  library(furrr)
  tryCatch(
    expr = plan(strategy = multisession(workers = threads)),
    error = function(x) {
    cat("something is wrong.\n")
    }
  )
  
  cor =
    rownames(microbiome_data) %>%
    furrr::future_map(function(name) {
      x = as.numeric(microbiome_data[name,])
      temp_partial_cor =
        purrr::map(
          as.data.frame(t(metabolome_data)),
          .f = function(y) {
            temp_data =
              data.frame(x = x,
                         y = y,
                         sample_info)
            temp_data$Gender[temp_data$Gender == 'F'] = 0
            temp_data$Gender[temp_data$Gender == 'M'] = 1
            temp_data$Gender = as.numeric(temp_data$Gender)
            
            temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
            temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
            temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
            temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
            temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
            
            temp_data$SSPG = as.numeric(temp_data$SSPG)
            temp_data$FPG = as.numeric(temp_data$FPG)
            
            ##only remain the subjects with less than 5 time points
            remain_subject_id =
              temp_data %>%
              dplyr::select(subject_id) %>%
              dplyr::group_by(subject_id) %>%
              dplyr::summarise(n = n()) %>%
              dplyr::filter(n >= 5) %>%
              pull(subject_id)
            
            temp_data =
              temp_data %>%
              dplyr::filter(subject_id %in% remain_subject_id)
            
            ##partial correlation
            if(method == "all"){
              cor_value1 =
                ppcor::pcor.test(
                  x = temp_data$x,
                  y = temp_data$y,
                  z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
                  method = "spearman"
                )
              
              result1 =
                c(cor = cor_value1$estimate,
                  p = cor_value1$p.value)
              
              if (is.na(result1[1])) {
                result1[1] = 0
              }
              
              if (is.na(result1[2])) {
                result1[2] = 1
              }
              
              cor_value2 =
                ppcor::pcor.test(
                  x = temp_data$x,
                  y = temp_data$y,
                  z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
                  method = "pearson"
                )
              
              result2 =
                c(cor = cor_value2$estimate,
                  p = cor_value2$p.value)
              
              if (is.na(result2[1])) {
                result2[1] = 0
              }
              
              if (is.na(result2[2])) {
                result2[2] = 1
              }
              
              list(result1 = result1,
                   result2 = result2)
              
            }else{
              cor_value =
                ppcor::pcor.test(
                  x = temp_data$x,
                  y = temp_data$y,
                  z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
                  method = method
                )
              
              result =
                c(cor = cor_value$estimate,
                  p = cor_value$p.value)
              
              if (is.na(result[1])) {
                result[1] = 0
              }
              
              if (is.na(result[2])) {
                result[2] = 1
              }
              list(result )
            }
          }
        ) 
        # do.call(rbind, .) %>%
        # as.data.frame()
      
      if(method == "all"){
        temp_partial_cor1 =
          temp_partial_cor %>% 
          purrr::map(function(x){
            x[[1]]
          }) %>% 
          do.call(rbind, .) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column(var = "metabolite") %>%
          dplyr::mutate(microbiome = name) %>%
          dplyr::select(microbiome, metabolite, dplyr::everything())
        
        temp_partial_cor1$p_adjust = p.adjust(temp_partial_cor1$p, method = "BH")
        
        temp_partial_cor2 =
          temp_partial_cor %>% 
          purrr::map(function(x){
            x[[2]]
          }) %>% 
          do.call(rbind, .) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column(var = "metabolite") %>%
          dplyr::mutate(microbiome = name) %>%
          dplyr::select(microbiome, metabolite, dplyr::everything())
        
        temp_partial_cor2$p_adjust = p.adjust(temp_partial_cor2$p, method = "BH")
        list(temp_partial_cor1 = temp_partial_cor1, 
             temp_partial_cor2 = temp_partial_cor2)
      }else{
        temp_partial_cor =
          temp_partial_cor %>%
          purrr::map(function(x){
            x[[1]]
          }) %>% 
          do.call(rbind, .) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column(var = "metabolite") %>%
          dplyr::mutate(microbiome = name) %>%
          dplyr::select(microbiome, metabolite, dplyr::everything())
        temp_partial_cor$p_adjust = p.adjust(temp_partial_cor$p, method = "BH")
        list(temp_partial_cor )
      }
    }, .progress = TRUE)
    # do.call(rbind, .) %>%
    # as.data.frame()
  
  if(method == "all"){
    cor1 = 
      cor %>% 
      purrr::map(function(x){
        x[[1]]
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    cor1$p_adjust2 =
      p.adjust(cor1$p, method = "BH")
    
    cor2 = 
      cor %>% 
      purrr::map(function(x){
        x[[2]]
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    cor2$p_adjust2 =
      p.adjust(cor2$p, method = "BH")
    list(cor1 = cor1, cor2 = cor2)    
  }else{
    cor = 
      cor %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    cor$p_adjust2 =
      p.adjust(cor$p, method = "BH")
    list(cor)
  }
  

}


####calculate the partial spearman correlation intra omics data
intra_partial_cor = function(expression_data,
                             sample_info) {
  library(future)
  library(furrr)
  plan(strategy = multisession(workers = 5))
  
  cor =
    1:(nrow(expression_data) - 1) %>%
    furrr::future_map(function(i) {
      # cat(i, " ")
      x = as.numeric(expression_data[i,])
      
      temp_cor =
        purrr::map(
          i:nrow(expression_data),
          .f = function(j) {
            y = as.numeric(expression_data[j,])
            temp_data =
              data.frame(x = x,
                         y = y,
                         sample_info)
            
            temp_data$Gender[temp_data$Gender == 'F'] = 0
            temp_data$Gender[temp_data$Gender == 'M'] = 1
            temp_data$Gender = as.numeric(temp_data$Gender)
            
            temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
            temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
            temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
            temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
            temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
            
            temp_data$SSPG = as.numeric(temp_data$SSPG)
            temp_data$FPG = as.numeric(temp_data$FPG)
            
            ##only remain the subjects with less than 5 time points
            remain_subject_id =
              temp_data %>%
              dplyr::select(subject_id) %>%
              dplyr::group_by(subject_id) %>%
              dplyr::summarise(n = n()) %>%
              dplyr::filter(n >= 5) %>%
              pull(subject_id)
            
            temp_data =
              temp_data %>%
              dplyr::filter(subject_id %in% remain_subject_id)
            
            ##partial correlation
            cor_value =
              ppcor::pcor.test(
                x = temp_data$x,
                y = temp_data$y,
                z = temp_data[, c("Gender", "Adj.age", "Ethnicity")],
                method = "spearman"
              )
            
            result =
              data.frame(
                name1 = rownames(expression_data)[i],
                name2 = rownames(expression_data)[2],
                cor = unname(cor_value$estimate),
                p = unname(cor_value$p.value)
              )
            
            if (is.na(result$cor)) {
              result$cor = 0
            }
            
            if (is.na(result$p)) {
              result$p = 1
            }
            result
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      
      temp_cor$p_adjust = p.adjust(temp_cor$p, method = "BH")
      temp_cor
    }, .progress = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  cor$p_adjust2 = p.adjust(cor$p, method = "BH")
  cor
}


####calculate the partial rmcorr between two omics data
rm_cor = function(microbiome_data,
                  metabolome_data,
                  sample_info) {
  cor =
    rownames(microbiome_data) %>%
    furrr::future_map(function(name) {
      cat(name, " ")
      x = as.numeric(microbiome_data[name,])
      temp_cor =
        purrr::map(
          as.data.frame(t(metabolome_data)),
          .f = function(y) {
            temp_data =
              data.frame(x = x,
                         y = y,
                         sample_info)
            temp_data$Gender[temp_data$Gender == 'F'] = 0
            temp_data$Gender[temp_data$Gender == 'M'] = 1
            temp_data$Gender = as.numeric(temp_data$Gender)
            
            temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
            temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
            temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
            temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
            temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
            
            temp_data$SSPG = as.numeric(temp_data$SSPG)
            temp_data$FPG = as.numeric(temp_data$FPG)
            
            ##only remain the subjects with less than 5 time points
            remain_subject_id =
              temp_data %>%
              dplyr::select(subject_id) %>%
              dplyr::group_by(subject_id) %>%
              dplyr::summarise(n = n()) %>%
              dplyr::filter(n >= 5) %>%
              pull(subject_id)
            
            temp_data =
              temp_data %>%
              dplyr::filter(subject_id %in% remain_subject_id)
            
            library(rmcorr)
            
            rm_cor =
              rmcorr(
                participant = "subject_id",
                measure1 = "x",
                measure2 = "y",
                dataset = temp_data
              )
            
            result =
              c(cor = rm_cor$r,
                p = rm_cor$p)
            
            if (is.na(result[1])) {
              result[1] = 0
            }
            
            if (is.na(result[2])) {
              result[2] = 1
            }
            result
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      temp_cor =
        temp_cor %>%
        tibble::rownames_to_column(var = "metabolite") %>%
        dplyr::mutate(microbiome = name) %>%
        dplyr::select(microbiome, metabolite, dplyr::everything())
      temp_cor$p_adjust = p.adjust(temp_cor$p, method = "BH")
      temp_cor
    }, .progress = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  cor$p_adjust2 =
    p.adjust(cor$p, method = "BH")
  cor
}




####calculate the rmcorr intra omics data
intra_rm_cor = function(expression_data,
                        sample_info) {
  library(future)
  library(furrr)
  plan(strategy = multisession(workers = 5))
  
  cor =
    1:(nrow(expression_data) - 1) %>%
    furrr::future_map(function(i) {
      x = as.numeric(expression_data[i, ])
      temp_cor =
        purrr::map(
          i:nrow(expression_data),
          .f = function(j) {
            y = as.numeric(expression_data[j, ])
            temp_data =
              data.frame(x = x,
                         y = y,
                         sample_info)
            
            temp_data$Gender[temp_data$Gender == 'F'] = 0
            temp_data$Gender[temp_data$Gender == 'M'] = 1
            temp_data$Gender = as.numeric(temp_data$Gender)
            
            temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
            temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
            temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
            temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
            temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
            
            temp_data$SSPG = as.numeric(temp_data$SSPG)
            temp_data$FPG = as.numeric(temp_data$FPG)
            
            ##only remain the subjects with less than 5 time points
            remain_subject_id =
              temp_data %>%
              dplyr::select(subject_id) %>%
              dplyr::group_by(subject_id) %>%
              dplyr::summarise(n = n()) %>%
              dplyr::filter(n >= 5) %>%
              pull(subject_id)
            
            temp_data =
              temp_data %>%
              dplyr::filter(subject_id %in% remain_subject_id)
            
            rm_cor =
              rmcorr(
                participant = "subject_id",
                measure1 = "x",
                measure2 = "y",
                dataset = temp_data
              )
            
            result =
              data.frame(
                name1 = rownames(expression_data)[i],
                name2 = rownames(expression_data)[2],
                cor = unname(rm_cor$r),
                p = rm_cor$p
              )
            
            if (is.na(result$cor)) {
              result$cor = 0
            }
            
            if (is.na(result$p)) {
              result$p = 1
            }
            result
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      
      temp_cor$p_adjust = p.adjust(temp_cor$p, method = "BH")
      temp_cor
    }, .progress = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  cor$p_adjust2 = p.adjust(cor$p, method = "BH")
  cor
}




individual_cor = 
  function(microbiome_data,
           metabolome_data,
           sample_info,
           method = c("spearman", "pearson")){
    method = match.arg(method)
    
    cor =
      sample_info$subject_id %>%
      unique() %>%
      purrr::map(function(temp_id) {
        idx = which(sample_info$subject_id == temp_id)
        temp_data1 = microbiome_data[,idx]
        temp_data2 = metabolome_data[,idx]

        cor_data =
        temp_data1 %>%
          t() %>%
          as.data.frame() %>%
        purrr::map(function(x){
          temp_cor =
          as.data.frame(t(temp_data2)) %>%
            purrr::map(function(y) {
              temp = data.frame(x, y)
              cor_test =
                cor.test(temp$x, temp$y, method = method)
              c(cor = unname(cor_test$estimate), p = cor_test$p.value)
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "metabolite")
          temp_cor$cor[is.na(temp_cor$cor)] = 0
          temp_cor$p[is.na(temp_cor$p)] = 1
          temp_cor
        })

        cor_data =
        purrr::map2(
          .x = cor_data,
          .y = names(cor_data),
          .f = function(x, y) {
            data.frame(x, microbiome = y) %>%
              dplyr::select(microbiome, metabolite, dplyr::everything())
          }
        ) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        cor_data
      })

    names(cor) =
      unique(sample_info$subject_id)

    ###calculate adjusted p values for each subject
    cor =
      purrr::map2(.x = names(cor),
                  .y = cor,
                  function(x, y) {
                    y$p_adjust = p.adjust(y$p, method = "BH")
                    y$subject_id = x
                    y
                  })

    cor =
      cor %>%
      do.call(rbind, .) %>%
      as.data.frame()

    rownames(cor) = NULL
    cor
  }









###three Transformation functions
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

arcsine <- function(p) {
  asin(sqrt(p))
}

logit <- function(y) {
  log(y / (1 - y))
}




modularity_plot = function(subnetworks){
  plot <- 
    ggplot(
      data.frame(index = 1:length(subnetworks$modularity),
                 modu = subnetworks$modularity, stringsAsFactors = FALSE),
      aes(index, modu) 
    ) +
    geom_vline(xintercept = which.max(subnetworks$modularity), 
               linetype = 2, colour = "#800000B2") + 
    labs(x = "Community analysis iteration", y = "Modularity") +
    geom_line(colour = "black") +
    # geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12))
  
  plot <-
    plot + 
    ggplot2::annotate(geom = "point", 
                      x = which.max(subnetworks$modularity),
                      y = max(subnetworks$modularity), 
                      size = 3, 
                      colour = "red") +
    annotate(geom = "text", 
             x = which.max(subnetworks$modularity),
             y = max(subnetworks$modularity), 
             label = paste("(",  which.max(subnetworks$modularity),
                           ",", 
                           max(subnetworks$modularity) %>% round(3),
                           ")"),
             size = 5,
             colour = "red"
    )
  
  plot
}


get_go_sim = 
  function(result, sim_cutoff = 0.7) {
  cc_sim_matrix <-
    simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "CC"],
                                      ont = "CC",
                                      measure = "Wang") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "name1") %>%
    tidyr::pivot_longer(cols = -name1,
                        names_to = "name2",
                        values_to = "sim") %>%
    dplyr::filter(name1 != name2) %>%
    dplyr::filter(sim > sim_cutoff)
  
  name <- apply(cc_sim_matrix, 1, function(x) {
    paste(sort(x[1:2]), collapse = "_")
  })
  
  cc_sim_matrix <-
    cc_sim_matrix %>%
    dplyr::mutate(name = name) %>%
    dplyr::arrange(name) %>%
    dplyr::distinct(name, .keep_all = TRUE) %>%
    dplyr::select(-name)
  
  bp_sim_matrix <-
    simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "BP"],
                                      ont = "BP",
                                      measure = "Wang") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "name1") %>%
    tidyr::pivot_longer(cols = -name1,
                        names_to = "name2",
                        values_to = "sim") %>%
    dplyr::filter(name1 != name2) %>%
    dplyr::filter(sim > sim_cutoff)
  
  name <- apply(bp_sim_matrix, 1, function(x) {
    paste(sort(x[1:2]), collapse = "_")
  })
  
  bp_sim_matrix <-
    bp_sim_matrix %>%
    dplyr::mutate(name = name) %>%
    dplyr::arrange(name) %>%
    dplyr::distinct(name, .keep_all = TRUE) %>%
    dplyr::select(-name)

  mf_sim_matrix <-
    simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "MF"],
                                      ont = "MF",
                                      measure = "Wang") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "name1") %>%
    tidyr::pivot_longer(cols = -name1,
                        names_to = "name2",
                        values_to = "sim") %>%
    dplyr::filter(name1 != name2) %>%
    dplyr::filter(sim > sim_cutoff)
  
  name <- apply(mf_sim_matrix, 1, function(x) {
    paste(sort(x[1:2]), collapse = "_")
  })
  
  mf_sim_matrix <-
    mf_sim_matrix %>%
    dplyr::mutate(name = name) %>%
    dplyr::arrange(name) %>%
    dplyr::distinct(name, .keep_all = TRUE) %>%
    dplyr::select(-name)
  
  sim_matrix = rbind(bp_sim_matrix,
                     cc_sim_matrix,
                     mf_sim_matrix)
  sim_matrix
}



get_go_cluster = function(result, sim_matrix, output_path = "."){
  edge_data <- 
    sim_matrix %>% 
    dplyr::rename(from = name1, to = name2)
  
  node_data <-
    result %>%
    as.data.frame() %>% 
    dplyr::select(ID, everything()) %>% 
    dplyr::rename(node = ID)
  
  graph <-
    tidygraph::tbl_graph(nodes = node_data,
                         edges = edge_data,
                         directed = FALSE) %>%
    dplyr::mutate(degree = tidygraph::centrality_degree())
  
  library(ggraph)
  library(igraph)
  library(tidygraph)
  
  subnetwork <-
    igraph::cluster_edge_betweenness(graph = graph,
                                     weights = abs(edge_attr(graph,
                                                             "sim")))
  cluster <-
    as.character(membership(subnetwork)) %>%
    purrr::map(function(x) {
      if (sum(x == as.character(membership(subnetwork))) == 1) {
        return("Other")
      } else{
        return(x)
      }
    }) %>%
    unlist()
  
  cluster1 <-
    purrr::map(cluster, function(x) {
      paste("Cluster", match(x, unique(cluster)[unique(cluster) != "Other"]))
    }) %>%
    unlist()
  
  cluster1[cluster1 == "Cluster NA"] <- "Other"
  
  graph <-
    graph %>%
    tidygraph::mutate(cluster = cluster1)
  
  result <-
    igraph::vertex_attr(graph) %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
    dplyr::arrange(ONTOLOGY, cluster, p.adjust)
  
  save(result, file = file.path(output_path, "result"))
  openxlsx::write.xlsx(
    result,
    file = file.path(output_path, "result.xlsx"),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  result_cluster <-
    result %>%
    dplyr::mutate(Count = as.numeric(Count)) %>% 
    plyr::dlply(.variables = .(cluster)) %>%
    purrr::map(function(x) {
      if (nrow(x) == 1) {
        x$cluster_annotation = x$Description
        return(x)
      }
      
      if (x$cluster[1] == 'Other') {
        x$cluster_annotation = x$Description
        return(x)
      }
      
      x$cluster_annotation = c(x$Description[which.min(as.numeric(x$p.adjust))],
                               x$Description[which.max(as.numeric(x$Count))]
                               ) %>% 
        unique() %>% 
        paste(., collapse = ";")
      
      x$node <-
        paste(x$node, collapse = ";")
      
      x$Description <-
        paste(x$Description, collapse = ";")
      
      x$BgRatio <-
        paste(x$BgRatio, collapse = ";")
      
      x$pvalue <- min(as.numeric(x$pvalue))
      x$p.adjust <- min(as.numeric(x$p.adjust))
      x$qvalue <- min(as.numeric(x$qvalue))
      x$geneID =
        x$geneID %>%
        stringr::str_split(pattern = "/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = '/')
      
      x$geneName =
        x$geneName %>%
        stringr::str_split(pattern = "/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = '/')
      
      x$Count <-
        sum(as.numeric(x$Count))
      
      x =
        x %>%
        dplyr::select(cluster, everything()) %>%
        dplyr::distinct(cluster, .keep_all = TRUE)
      
      x
      
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(cluster_annotation, everything())
  
  save(result_cluster, file = file.path(output_path, "result_cluster"))
  openxlsx::write.xlsx(
    result_cluster,
    file = file.path(output_path, "result_cluster.xlsx"),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  ####clustered different GO terms
  for(ont in c("CC", "MF", "BP")) {
    cat(ont, " ")
    plot <-
      graph %>%
      tidygraph::filter(ONTOLOGY == ont) %>%
      tidygraph::filter(cluster != "Other") %>%
      ggraph(layout = 'fr',
             circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        strength = 1,
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = cluster,
            size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(cluster == "Other", NA, Description)
      ),
      size = 3,
      repel = TRUE) +
      guides(fill = guide_legend(ncol = 1)) +
      # scale_fill_manual(values = c(
      #   "Other" = "grey"
      # )) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    ggsave(
      plot,
      filename = file.path(output_path,
                           paste(ont, "sim_plot.pdf", sep = "_")),
      width = 9,
      height = 7
    )
  }
  
  for(temp_cluster in unique(cluster1)) {
    cat(temp_cluster, " ")
    plot1 <-
      graph %>%
      tidygraph::filter(cluster == temp_cluster) %>%
      ggraph(layout = 'fr',
             circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        strength = 1,
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = cluster,
            size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(cluster == "Other", NA, Description)
      ),
      size = 3,
      repel = TRUE) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "left",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    plot2 <-
      result %>%
      dplyr::filter(cluster == temp_cluster) %>%
      dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(Description = factor(Description, levels = Description)) %>%
      ggplot(aes(p.adjust, Description)) +
      geom_bar(stat = "identity") +
      geom_text(aes(x = 0, Description, label = Description),
                hjust = 0,
                size = 5) +
      theme_bw() +
      labs(y = "", x = "-log10(FDR adjusted P value)") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    library(patchwork)
    plot =
      plot1 + plot2 + patchwork::plot_layout(nrow = 1)
    
    ggsave(
      plot,
      filename = file.path(output_path,
                           paste(temp_cluster, "sim_plot.pdf", sep = "_")),
      width = 14,
      height = 7
    )
  }
}


output_go_structure = function(result_cluster, output_path = ".") {
  unique(result_cluster$cluster) %>%
    purrr::map(
      .f = function(x) {
        cat(x, " ")
        if (x == "Other") {
          return(NULL)
        }
        
        temp_id =
          result_cluster %>%
          dplyr::filter(cluster == x) %>%
          dplyr::pull(node) %>%
          stringr::str_split(";") %>%
          `[[`(1) %>%
          pRoloc::goIdToTerm(keepNA = FALSE) %>%
          data.frame(id = ., class = "YES") %>%
          tibble::rownames_to_column(var = "name")
        
        temp_plot =
          GOSim::getGOGraph(term = temp_id$name, prune = Inf) %>%
          igraph.from.graphNEL() %>%
          tidygraph::as_tbl_graph() %>%
          tidygraph::left_join(temp_id, by = "name") %>%
          dplyr::mutate(class = case_when(is.na(class) ~ "NO",
                                          TRUE ~ class))
        
        plot =
          temp_plot %>%
          ggraph(layout = 'kk',
                 circular = FALSE) +
          geom_edge_link(
            color = ggsci::pal_aaas()(n = 10)[1],
            alpha = 1,
            arrow = grid::arrow(
              angle = 10,
              length = unit(0.2, "inches"),
              type = "closed"
            ),
            show.legend = FALSE
          ) +
          geom_node_point(
            aes(fill = class),
            shape = 21,
            alpha = 1,
            size = 6,
            show.legend = FALSE
          ) +
          geom_node_text(aes(
            x = x,
            y = y,
            label = ifelse(class == "YES", id, NA)
          ),
          size = 3,
          repel = TRUE) +
          scale_fill_manual(values = c('YES' = "red", 'NO' = "white")) +
          ggraph::theme_graph() +
          theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "left",
            legend.background = element_rect(fill = "transparent", color = NA)
          )
        
        ggsave(
          plot,
          filename = file.path(output_path,
                               paste(x, "_GO graph.pdf", sep = "")),
          width = 7,
          height = 7
        )
      }
    )
}






#' sparcc wrapper
#'
#' A reimplementation of SparCC algorithm (Friedman et Alm 2012, PLoS Comp Bio, 2012).
#' @param data Community count data matrix
#' @param iter Number of iterations in the outer loop
#' @param inner_iter Number of iterations in the inner loop
#' @param th absolute value of correlations below this threshold are considered zero by the inner SparCC loop.
#' @seealso \code{\link{sparccboot}}
#' @export
sparcc <- function(data, iter=20, inner_iter=10, th=.1) {
  ##
  #  without all the 'frills'
  sparccs <- lapply(1:iter, function(i)
    sparccinner(t(apply(data, 1, norm_diric)),
                iter=inner_iter, th=th))
  # collect
  cors <- array(unlist(lapply(sparccs, function(x) x$Cor)),
                c(ncol(data),ncol(data),iter))
  corMed <- apply(cors, 1:2, median)
  covs <- array(unlist(lapply(sparccs, function(x) x$Cov)),
                c(ncol(data),ncol(data),iter))
  covMed <- apply(covs, 1:2, median)
  covMed <- cor2cov(corMed, sqrt(diag(covMed)))
  list(Cov=covMed, Cor=corMed)
}

#' Bootstrap SparCC
#'
#' Get bootstrapped estimates of SparCC correlation coefficients. To get empirical p-values, pass this output to \code{pval.sparccboot}.
#'
#' @param data Community count data
#' @param sparcc.params named list of parameters to pass to \code{sparcc}
#' @param statisticboot function which takes data and bootstrap sample indices and results the upper triangle of the bootstapped correlation matrix
#' @param statisticperm function which takes data and permutated sample indices and results the upper triangle of the null correlation matrix
#' @param R number of bootstraps
#' @param ncpus number of cores to use for parallelization
#' @param ... additional arguments that are passed to \code{boot::boot}
#' @export
sparccboot <- function(data, sparcc.params=list(),
                       statisticboot=function(data, indices) triu(do.call("sparcc",
                                                                          c(list(data[indices,,drop=FALSE]), sparcc.params))$Cor),
                       statisticperm=function(data, indices) triu(do.call("sparcc",  c(list(apply(data[indices,], 2, sample)), sparcc.params))$Cor),
                       R, ncpus=1, ...) {
  
  if (!requireNamespace('boot', quietly=TRUE))
    stop('\'boot\' package is not installed')
  
  res     <- boot::boot(data, statisticboot, R=R, parallel="multicore", ncpus=ncpus, ...)
  null_av <- boot::boot(data, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=ncpus)
  class(res) <- 'list'
  structure(c(res, list(null_av=null_av)), class='sparccboot')
}

#' SparCC p-vals
#'
#' Get empirical p-values from bootstrap SparCC output.
#'
#' @param x output from \code{sparccboot}
#' @param sided type of p-value to compute. Only two sided (sided="both") is implemented.
#' @export
pval.sparccboot <- function(x, sided='both') {
  # calculate 1 or 2 way pseudo p-val from boot object
  # Args: a boot object
  if (sided != "both") stop("only two-sided currently supported")
  nparams  <- ncol(x$t)
  tmeans   <- colMeans(x$null_av$t)
  #    check to see whether correlations are unstable -- confirm
  #    that sample correlations are in 95% confidence interval of
  #    bootstrapped samples
  niters   <- nrow(x$t)
  ind95    <- max(1,round(.025*niters)):round(.975*niters)
  boot_ord <- apply(x$t, 2, sort)
  boot_ord95 <- boot_ord[ind95,]
  outofrange <- unlist(lapply(1:length(x$t0), function(i) {
    aitvar <- x$t0[i]
    range  <- range(boot_ord95[,i])
    range[1] > aitvar || range[2] < aitvar
  }))
  # calc whether center of mass is above or below the mean
  bs_above <- unlist(lapply(1:nparams, function(i)
    length(which(x$t[, i] > tmeans[i]))))
  is_above <- bs_above > x$R/2
  cors <- x$t0
  #    signedAV[is_above] <- -signedAV[is_above]
  pvals    <- ifelse(is_above, 2*(1-bs_above/x$R), 2*bs_above/x$R)
  pvals[pvals > 1]  <- 1
  pvals[outofrange] <- NaN
  list(cors=cors, pvals=pvals)
}



#' @noRd
sparccinner <- function(data.f, T=NULL, iter=10, th=0.1) {
  if (is.null(T))   T  <- av(data.f)
  res.bv <- basis_var(T)
  Vbase  <- res.bv$Vbase
  M      <- res.bv$M
  cbase  <- C_from_V(T, Vbase)
  Cov    <- cbase$Cov
  Cor    <- cbase$Cor
  
  ## do iterations here
  excluded <- NULL
  for (i in 1:iter) {
    res.excl <- exclude_pairs(Cor, M, th, excluded)
    M <- res.excl$M
    excluded <- res.excl$excluded
    if (res.excl$break_flag) break
    res.bv <- basis_var(T, M=M, excluded=excluded)
    Vbase  <- res.bv$Vbase
    M      <- res.bv$M
    K <- M
    diag(K) <- 1
    cbase  <- C_from_V(T, Vbase)
    Cov    <- cbase$Cov
    Cor    <- cbase$Cor
  }
  list(Cov=Cov, Cor=Cor, i=i, M=M, excluded=excluded)
}

#' @noRd
exclude_pairs <- function(Cor, M, th=0.1, excluded=NULL) {
  # exclude pairs with high correlations
  break_flag <- FALSE
  C_temp <- abs(Cor - diag(diag(Cor)) )  # abs value / remove diagonal
  if (!is.null(excluded)) C_temp[excluded] <- 0 # set previously excluded correlations to 0
  exclude <- which(abs(C_temp - max(C_temp)) < .Machine$double.eps*100)[1:2]
  if (max(C_temp) > th)  {
    i <- na.exclude(arrayInd(exclude, c(nrow(M), ncol(M)))[,1])
    M[i,i] <- M[i,i] - 1
    excluded_new <- c(excluded, exclude)
  } else {
    excluded_new <- excluded
    break_flag   <- TRUE
  }
  list(M=M, excluded=excluded_new, break_flag=break_flag)
}

#' @noRd
basis_cov <- function(data.f) {
  # data.f -> relative abundance data
  # OTUs in columns, samples in rows (yes, I know this is transpose of normal)
  # first compute aitchison variation
  T <- av(data.f)
  res.bv <- basis_var(T)
  Vbase  <- res.bv$Vbase
  M      <- res.bv$M
  cbase  <- C_from_V(T, Vbase)
  Cov    <- cbase$Cov
  Cor    <- cbase$Cor
  list(Cov=Cov, M=M)
}

#' @noRd
basis_var <- function(T, CovMat = matrix(0, nrow(T), ncol(T)),
                      M = matrix(1, nrow(T), ncol(T)) + (diag(ncol(T))*(ncol(T)-2)),
                      excluded = NULL, Vmin=1e-4) {
  
  if (!is.null(excluded)) {
    T[excluded] <- 0
    #   CovMat[excluded] <- 0
  }
  Ti     <- matrix(rowSums(T))
  CovVec <- matrix(rowSums(CovMat - diag(diag(CovMat)))) # row sum of off diagonals
  M.I <- tryCatch(solve(M), error=function(e) MASS::ginv(M))
  Vbase <- M.I %*% (Ti + 2*CovVec)
  Vbase[Vbase < Vmin] <- Vmin
  list(Vbase=Vbase, M=M)
}

C_from_V <- function(T, Vbase) {
  J      <- matrix(1, nrow(T), ncol(T))
  Vdiag  <- diag(c(Vbase))
  CovMat <- .5*((J %*% Vdiag) + (Vdiag %*% J) - T)
  CovMat <- (CovMat + t(CovMat))/2  # enforce symmetry
  # check that correlations are within -1,1
  CorMat <- cov2cor(CovMat)
  CorMat[abs(CorMat) > 1] <- sign(CorMat[abs(CorMat) > 1])
  CovMat <- cor2cov(CorMat, sqrt(as.vector(Vbase)))
  list(Cov=CovMat, Cor=CorMat)
}

#' @noRd
av <- function(data) {
  cov.clr <- cov(clr(data))
  J <- matrix(1, ncol(data), ncol(data))
  (J %*% diag(diag(cov.clr))) + (diag(diag(cov.clr)) %*% J) - (2*cov.clr)
}


#' @importFrom VGAM rdiric
#' @noRd
norm_diric   <- function(x, rep=1) {
  dmat <- VGAM::rdiric(rep, x+1)
  norm_to_total(colMeans(dmat))
}


norm_to_total <- function(x) x/sum(x)


#' Convert a symmetric correlation matrix to a covariance matrix
#' given the standard deviation
#'
#' @param cor a symmetric correlation matrix
#' @param sds standard deviations of the resulting covariance.
#' @return Covariance matrix of sample dimension as cor
#' @export
cor2cov <- function(cor, sds) {
  if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
  cor * sds * rep(sds, each=nrow(cor))
}