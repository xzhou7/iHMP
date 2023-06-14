
library(plyr)
library(dplyr)

source("./common/merge_metadata.r")
source("./common/match_datasets.r")
source("./common/theme_nature.r")

source("./common/load_bugs.r")
source("./common/load_pathways.r")
source("./common/load_metabolites.r")
source("./common/load_proteins.r")
source("./common/load_fecalcal.r")



subway_venn <- function(diagnoses, lim=1000) {
    df <- data.frame()

    add_bar <- function(name, pcl_list) {
        pcl_list <- lapply(pcl_list, pcl.filter.s, diagnosis %in% diagnoses, enclos = parent.frame())
        strict    <- nrow(match_datasets(pcl_list, lenience=0, matching=T))
        lenience2 <- nrow(match_datasets(pcl_list, lenience=2, matching=T))
        lenience4 <- nrow(match_datasets(pcl_list, lenience=4, matching=T))
        df <<- rbind(df,
                     data.frame(grp = name, type = "strict",    count = strict),
                     data.frame(grp = name, type = "lenience2", count = lenience2-strict),
                     data.frame(grp = name, type = "lenience4", count = lenience4-lenience2)
        )
    }

    add_bar("MGX+MTX", list(pwy.dna.unstrat.pcl, pwy.rna.unstrat.pcl))
    add_bar("MGX+FC", list(pwy.dna.unstrat.pcl, fecalcal.pcl))
    add_bar("MGX+MTX+FC", list(pwy.dna.unstrat.pcl, pwy.rna.unstrat.pcl, fecalcal.pcl))
    add_bar("MGX+MTX+VX", list(pwy.dna.unstrat.pcl, pwy.rna.unstrat.pcl, viruses.pcl))
    add_bar("MGX+MPX", list(pwy.dna.unstrat.pcl, proteins.pcl))
    add_bar("MGX+MBX", list(pwy.dna.unstrat.pcl, metabolites.pcl))
    add_bar("MGX+MPX+MBX", list(pwy.dna.unstrat.pcl, proteins.pcl, metabolites.pcl))
    add_bar("MGX+MTX+MPX", list(pwy.dna.unstrat.pcl, pwy.rna.unstrat.pcl, proteins.pcl))
    add_bar("MGX+MTX+MBX", list(pwy.dna.unstrat.pcl, pwy.rna.unstrat.pcl, metabolites.pcl))
    add_bar("ALL", list(pwy.dna.unstrat.pcl, pwy.rna.unstrat.pcl, proteins.pcl, metabolites.pcl, fecalcal.pcl))

    df$grp.f <- factor(df$grp, levels=rev(df$grp[!duplicated(df$grp)]))
    df$type <- factor(df$type, levels=c("lenience4", "lenience2", "strict"))

    library(ggplot2)
    library(cowplot)
    ggp <- ggplot(data=df, aes(x=grp.f, y=count)) +
        geom_bar(aes(fill=type), color="black", stat="identity", width=0.8, lwd=0.1) +
        scale_y_continuous(position = "top", limits=c(0, lim), expand=c(0, 0)) +
        coord_flip() + theme_nature() +
        scale_fill_manual(values=c(strict="#7ba8cc", lenience2="#e7e09f", lenience4="#e0848e"),
                          breaks=c("strict", "lenience2", "lenience4"),
                          labels=c("Strict", "± 2 weeks", "± 4 weeks")) +
        guides(fill=guide_legend(title="Strictness")) +
        stat_summary(fun.data = function(x){
            return(c(y = sum(x) + 20, label = sum(x)))
        }, geom = "text", position = position_dodge(width = 0.75), hjust="left", vjust="middle", size=1.6) +
        stat_summary(fun.data = function(x){
            return(c(y = x[1] - 20, label = x[1]))
        }, geom = "text", position = position_dodge(width = 0.75), hjust="right", vjust="middle", color="white", size=1.6) +
        ylab("Timepoints") + xlab(NULL)

    return (ggp)
}


ggp <- subway_venn(c("CD", "UC", "nonIBD"))
pdf("./overview/sample_size_venn_barplots.pdf", 2.1, 1.54)
print(ggp)
dev.off()


library(egg)

ggp <- ggarrange(
    subway_venn("CD", 500),
    subway_venn("UC", 500),
    subway_venn("nonIBD", 500),
    ncol=1)
pdf("./overview/sample_size_venn_barplots_bydisease.pdf", 2.359, 4.839)
print(ggp)
dev.off()
