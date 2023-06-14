
library(ggplot2)
source("./common/theme_nature.r")


source("./common/load_bugs.r")
source("./common/load_kos.r")
source("./common/load_metabolites.r")
source("./common/load_proteins.r")
source("./common/load_biopsies.r")
source("./common/load_diet.r")

source("./common/match_datasets.r")

doa_aod_plot <- function(name1, pcl1, method1, name2, pcl2, method2, lenience=2) {
    pcl.matched <- match_datasets(list(pcl1, pcl2), lenience = lenience)

    library(vegan)
    D1 <- vegdist(pcl.matched[[1]]$x, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(pcl.matched[[2]]$x, method=method2, binary=method2=="jaccard")

    unqSubj <- as.character(unique(pcl.matched[[1]]$meta$subject))
    subject <- factor(pcl.matched[[1]]$meta$subject, levels=unqSubj)

    # AOD
    library(ade4)
    Daod1 <- matrix(0, length(unqSubj), length(unqSubj))
    rownames(Daod1) <- unqSubj
    colnames(Daod1) <- unqSubj
    Daod2 <- Daod1
    D1m <- as.matrix(D1)
    D2m <- as.matrix(D2)
    for (i in seq_along(unqSubj)) {
        for (j in seq_along(unqSubj)) {
            if (i > j) {
                mu1 <- mean(D1m[subject==unqSubj[i], subject==unqSubj[j]])
                Daod1[i,j] <- mu1
                Daod1[j,i] <- mu1
                mu2 <- mean(D2m[subject==unqSubj[i], subject==unqSubj[j]])
                Daod2[i,j] <- mu2
                Daod2[j,i] <- mu2
            }
        }
    }
    Daod1 <- as.dist(Daod1)
    Daod2 <- as.dist(Daod2)
    Daod1m <- as.matrix(Daod1)
    Daod2m <- as.matrix(Daod2)

    # DOA
    library(ade4)
    ave1.pcl <- pcl.group.s(pcl.matched[[1]], subj, avg=T)
    ave2.pcl <- pcl.group.s(pcl.matched[[2]], subj, avg=T)

    library(vegan)
    Ddoa1 <- vegdist(ave1.pcl$x, method=method1, binary=method1=="jaccard")
    Ddoa2 <- vegdist(ave2.pcl$x, method=method2, binary=method2=="jaccard")
    Ddoa1m <- as.matrix(Ddoa1)
    Ddoa2m <- as.matrix(Ddoa2)

    pdf(sprintf("./overview/mantel_doa_aod_debug_%s_%s.pdf", name1, name2), 9, 8)
    for (i in 1:100) {
        subjs <- sample(unqSubj, 3)
        keepS <- subject %in% subjs
        keepSd <- unclass(as.dist(outer(keepS, keepS, FUN="&")))>0
        same_subj <- unclass(as.dist(outer(subject, subject, FUN="==")))>0
        df_ds <- data.frame(
            D1 = unclass(D1)[keepSd & !same_subj],
            D2 = unclass(D2)[keepSd & !same_subj])
        df_l <- data.frame()
        for (s1i in seq_along(subjs)) {
            for (s2i in seq_along(subjs)) {
                if (s1i > s2i) {
                    s1 <- subjs[s1i]
                    s2 <- subjs[s2i]
                    iss12 <- (subject == s1) | (subject == s2)
                    keepSd12 <- unclass(as.dist(outer(iss12, iss12, FUN="&")))>0
                    p1 <- unclass(D1)[keepSd12 & !same_subj]
                    p2 <- unclass(D2)[keepSd12 & !same_subj]
                    df_l <- rbind(df_l,
                        data.frame(
                            D1 = c(rep(Ddoa1m[unqSubj == s1, unqSubj == s2], length(p1)), p1),
                            D2 = c(rep(Ddoa2m[unqSubj == s1, unqSubj == s2], length(p2)), p2),
                            pt = rep(seq_along(p1)+nrow(df_l), 2),
                            doa_aod = rep("DOA", 2*length(p1))),
                        data.frame(
                            D1 = c(rep(Daod1m[unqSubj == s1, unqSubj == s2], length(p1)), p1),
                            D2 = c(rep(Daod2m[unqSubj == s1, unqSubj == s2], length(p2)), p2),
                            pt = rep(seq_along(p1)+nrow(df_l)+length(p1), 2),
                            doa_aod = rep("AOD", 2*length(p2))))
                }
            }
        }

        ggp <- ggplot() +
            geom_line(data=df_l, aes(x=D1, y=D2, group=pt, color=doa_aod)) +
            scale_color_manual(values=c(AOD="darkorange2", DOA="deepskyblue2")) +
            geom_point(data=df_ds, aes(x=D1, y=D2)) +
            xlab(name1) + ylab(name2)

        print(ggp)
    }
    dev.off()
}


interindividual_mantel_test_aod <- function(name1, pcl1, method1, name2, pcl2, method2, Nperms=999, lenience = 2, Nbs=Nperms, ...) {
    pcl.matched <- match_datasets(list(pcl1, pcl2), lenience = lenience)

    library(vegan)
    D1 <- vegdist(pcl.matched[[1]]$x, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(pcl.matched[[2]]$x, method=method2, binary=method2=="jaccard")

    # Build matrices of the averages-of-dissimilarities
    library(ade4)
    subject <- pcl.matched[[1]]$meta$subject
    unqSubj <- unique(subject)
    Da1 <- matrix(0, length(unqSubj), length(unqSubj))
    rownames(Da1) <- unqSubj
    colnames(Da1) <- unqSubj
    Da2 <- Da1
    D1 <- as.matrix(D1)
    D2 <- as.matrix(D2)
    for (i in seq_along(unqSubj)) {
        for (j in seq_along(unqSubj)) {
            if (i > j) {
                mu1 <- mean(D1[subject==unqSubj[i], subject==unqSubj[j]])
                Da1[i,j] <- mu1
                Da1[j,i] <- mu1
                mu2 <- mean(D2[subject==unqSubj[i], subject==unqSubj[j]])
                Da2[i,j] <- mu2
                Da2[j,i] <- mu2
            }
        }
    }
    Da1 <- as.dist(Da1)
    Da2 <- as.dist(Da2)

    mantelbootstrap <- function(m1, m2, nrepet) {
        s1 <- unclass(m1)
        s2 <- unclass(m2)
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        bss <- apply(permi, 1, function(i) {
            bsi <- sample(seq_along(s1), length(s1), replace=T)
            return (cor(s1[bsi], s2[bsi]))
        })
        return (bss)
    }

    mt <- mantel.rtest(Da1, Da2, nrepet=Nperms)
    mt$bootstraps <- mantelbootstrap(Da1, Da2, nrepet=Nbs)
    return (mt)
}

interindividual_mantel_test_doa <- function(name1, pcl1, method1, name2, pcl2, method2, Nperms=999, lenience = lenience, Nbs=Nperms, ...) {
    pcl.matched <- match_datasets(list(pcl1, pcl2), lenience = lenience)

    # Build matrices of the dissimilarities-of-averages
    library(ade4)
    unqSubj <- unique(pcl.matched[[1]]$meta$subject)
    subj <- factor(pcl.matched[[1]]$meta$subject, levels=unqSubj)
    ave1.pcl <- pcl.group.s(pcl.matched[[1]], subj, avg=T)
    ave2.pcl <- pcl.group.s(pcl.matched[[2]], subj, avg=T)

    library(vegan)
    D1 <- vegdist(ave1.pcl$x, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(ave2.pcl$x, method=method2, binary=method2=="jaccard")

    mantelbootstrap <- function(m1, m2, nrepet) {
        s1 <- unclass(m1)
        s2 <- unclass(m2)
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        bss <- apply(permi, 1, function(i) {
            bsi <- sample(seq_along(s1), length(s1), replace=T)
            return (cor(s1[bsi], s2[bsi]))
        })
        return (bss)
    }

    mt <- mantel.rtest(D1, D2, nrepet=Nperms)
    mt$bootstraps <- mantelbootstrap(D1, D2, nrepet=Nbs)
    return (mt)
}

intraindividual_mantel_test <- function(name1, pcl1, method1, name2, pcl2, method2, lenience, Nperms=999, Nbs=Nperms, ...) {
    pcl.matched <- match_datasets(list(pcl1, pcl2), lenience = lenience)

    library(vegan)
    D1 <- vegdist(pcl.matched[[1]]$x, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(pcl.matched[[2]]$x, method=method2, binary=method2=="jaccard")

    # Code based on ade4:::mantel.rtest

    library(permute)
    library(ade4)
    permhow <- how(blocks=pcl.matched[[1]]$meta$subject)
    intrasubjectD <- unclass(as.dist(outer(
        as.character(pcl.matched[[1]]$meta$subject),
        as.character(pcl.matched[[1]]$meta$subject), FUN="=="))) > 0

    if (sum(intrasubjectD) < 100) {
        return (list(obs = NA, pvalue = NA, bootstraps = NA))
    }
    permutedist <- function(m, i) {
        w0 <- permute(i, pcl.matched[[1]]$ns, permhow)
        m <- as.matrix(m)
        return(as.dist(m[w0, w0]))
    }
    mantelnoneuclid <- function(m1, m2, nrepet) {
        obs <- cor(unclass(m1)[intrasubjectD], unclass(m2)[intrasubjectD])
        if (nrepet == 0)
            return(obs)
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        perm <- apply(permi, 1, function(i)
            cor(unclass(m1)[intrasubjectD],
                unclass(permutedist(m2, i))[intrasubjectD]))
        w <- as.randtest(obs = obs, sim = perm, call = match.call(),
                         subclass = "mantelrtest", ...)
        return(w)
    }

    mantelbootstrap <- function(m1, m2, nrepet) {
        s1 <- unclass(m1)[intrasubjectD]
        s2 <- unclass(m2)[intrasubjectD]
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        bss <- apply(permi, 1, function(i) {
            bsi <- sample(seq_along(s1), length(s1), replace=T)
            return (cor(s1[bsi], s2[bsi]))
        })
        return (bss)
    }

    mt <- mantelnoneuclid(D1, D2, nrepet=Nperms)
    mt$bootstraps <- mantelbootstrap(D1, D2, nrepet=Nbs)
    return (mt)
}

pcls <- list(
    "Taxonomy"    = bugs.pcl %>% pcl.only(rank="s") %>% pcl.normalize,
    "KOs (DNA)"   = ko.dna.unstrat.pcl %>% pcl.normalize,
    "KOs (RNA)"   = ko.rna.unstrat.pcl %>% pcl.normalize,
    "KOs (Protein)"=proteins.kos.pcl %>% pcl.normalize,
    "Metabolites" = metabolites.pcl.nrm,
    #"Viromics"    = viruses.pcl %>% pcl.filter.s(any(x>0)),
    "Biopsy 16S"  = biopsy_16s.pcl %>% pcl.normalize,
    "Biopsy HTX"  = biopsy_htx.counts.pcl %>% pcl.normalize,
    "Diet"        = hmp2_diet_num.pcl %>% pcl.filter.s(!any(is.na(x)))
)
distance_method <- c(
    "Taxonomy"    = "bray",
    "KOs (DNA)"   = "bray",
    "KOs (RNA)"   = "bray",
    "KOs (Protein)"="bray",
    "Metabolites" = "bray",
    "Viromics"    = "jaccard",
    "Biopsy 16S"  = "bray",
    "Biopsy HTX"  = "bray",
    "Diet"        = "manhattan"
)
match_lenience <- c(
    "Taxonomy"    = 0,
    "KOs (DNA)"   = 0,
    "KOs (RNA)"   = 2,
    "KOs (Protein)"=2,
    "Metabolites" = 2,
    "Viromics"    = 2,
    "Biopsy 16S"  = 4,
    "Biopsy HTX"  = 4,
    "Diet"        = 0
)

# pb <- txtProgressBar(min = 0, max = length(pcls)*length(pcls), style = 3)
# for (i in seq_along(pcls)) {
#     for (j in seq_along(pcls)) {
#         if (i > j) {
#             lenience <- max(match_lenience[names(pcls)[i]], match_lenience[names(pcls)[j]])
#
#             doa_aod_plot(
#                 names(pcls)[i], pcls[[i]], distance_method[names(pcls)[i]],
#                 names(pcls)[j], pcls[[j]], distance_method[names(pcls)[j]],
#                 lenience=lenience)
#         }
#
#         setTxtProgressBar(pb, (i-1)*length(pcls) + j)
#     }
# }
# close(pb)


M <- matrix(NA, nrow=length(pcls), ncol=length(pcls))
rownames(M) <- names(pcls)
colnames(M) <- names(pcls)

mt_intra <- list(C=M, Cil=M, Ciu=M, P=M)
mt_inter_aod <- mt_intra
mt_inter_doa <- mt_intra
Nperms <- 4999

pb <- txtProgressBar(min = 0, max = length(pcls)*length(pcls), style = 3)
for (i in seq_along(pcls)) {
    for (j in seq_along(pcls)) {
        if (i > j && is.na(mt_inter_doa$C[i,j])) {
            lenience <- max(match_lenience[names(pcls)[i]], match_lenience[names(pcls)[j]])

            mt <- intraindividual_mantel_test(
                names(pcls)[i], pcls[[i]], distance_method[names(pcls)[i]],
                names(pcls)[j], pcls[[j]], distance_method[names(pcls)[j]],
                Nperms=Nperms, lenience=lenience)
            mt_intra$C[i,j] <- mt$obs
            mt_intra$Cil[i,j] <- quantile(mt$bootstraps, 0.025, na.rm=T)
            mt_intra$Ciu[i,j] <- quantile(mt$bootstraps, 0.975, na.rm=T)
            mt_intra$P[i,j] <- mt$pvalue

            # mt <- interindividual_mantel_test_aod(
            #     names(pcls)[i], pcls[[i]], distance_method[names(pcls)[i]],
            #     names(pcls)[j], pcls[[j]], distance_method[names(pcls)[j]],
            #     Nperms=Nperms, lenience=lenience)
            # mt_inter_aod$C[i,j] <- mt$obs
            # mt_inter_aod$Cil[i,j] <- quantile(mt$bootstraps, 0.025)
            # mt_inter_aod$Ciu[i,j] <- quantile(mt$bootstraps, 0.975)
            # mt_inter_aod$P[i,j] <- mt$pvalue

            mt <- interindividual_mantel_test_doa(
                names(pcls)[i], pcls[[i]], distance_method[names(pcls)[i]],
                names(pcls)[j], pcls[[j]], distance_method[names(pcls)[j]],
                Nperms=Nperms, lenience=lenience)
            mt_inter_doa$C[i,j] <- mt$obs
            mt_inter_doa$Cil[i,j] <- quantile(mt$bootstraps, 0.025, na.rm=T)
            mt_inter_doa$Ciu[i,j] <- quantile(mt$bootstraps, 0.975, na.rm=T)
            mt_inter_doa$P[i,j] <- mt$pvalue
        }

        setTxtProgressBar(pb, (i-1)*length(pcls) + j)
    }
}
close(pb)

save(mt_intra, file="./overview/mantel_test_intra.RData")
#save(mt_inter_aod, file="./overview/mantel_test_inter_aod.RData")
save(mt_inter_doa, file="./overview/mantel_test_inter_doa.RData")

#load(file="./overview/mantel_test_intra.RData")
#load(file="./overview/mantel_test_inter_doa.RData")

# Plots
library(RColorBrewer)
manteltest_plot <- function(O, P, Ocil, Ociu, title, stars=F) {

    library(reshape2)
    library(viridis)
    df <- melt(O)
    colnames(df) <- c("V1", "V2", "obs")
    df$V2 <- factor(df$V2, levels=rev(levels(df$V2)))
    df$obs <- df$obs^2
    if (!missing(P)) {
        Padj <- matrix(p.adjust(P, method="fdr"), nrow=nrow(P))
        df$padj <- melt(Padj)$value
    }
    if (!missing(Ocil)) {
        df$cil <- melt(Ocil)$value^2
        df$cil[melt(Ocil)$value < 0] <- 0
        df$ciu <- melt(Ociu)$value^2
    }

    # Increase the contrast of the color scale where it matters
    alpha <- 1.9
    beta <- 0.1
    colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
    df$lblcolor <- ifelse(qbeta(df$obs, alpha, beta) < 0.8, "black", "white")

    ggp <- ggplot(data=df, aes(x=V1, y=V2)) +
        geom_tile(aes(fill=obs)) +
        scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                             values=colorvalues,
                             na.value="white", limits=c(0, 1), name="Variance Explained")
    #ggp
    nudge_y <- 0
    if (!missing(P)) {
        ggp <- ggp + if (stars) {
            nudge_y <- -0.14
            geom_text(aes(label=ifelse(padj<=0.001, "***", ifelse(padj<=0.01, "**", ifelse(padj<=0.05, "*", ""))),
                          color=lblcolor),
                      size=2, nudge_y=0.15, fontface="bold")
        } else {
            geom_text(aes(label=ifelse(is.na(padj), "", sprintf("FDR p\n%.2g", padj)),
                          color=lblcolor), size=1.7)
        }
    }
    ggp <- ggp + if (!missing(Ocil)) {
        geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%\n[%.1f%% - %.1f%%]", 100*obs, 100*cil, 100*ciu)),
                      color=lblcolor),
                  size=1.7, nudge_y=nudge_y)
    } else {
        geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%", 100*obs)),
                      color=lblcolor),
                  size=1.7, nudge_y=nudge_y)
    }
    ggp <- ggp +
        geom_text(aes(label=ifelse(V1!=V2, ifelse(is.na(obs), "N/A", ""), "")), size=1.7, color="gray") +
        theme_nature() +
        theme(axis.text.x=element_text(angle=-17, hjust=0.05),
              legend.position = "left", axis.ticks.y = element_blank()) +
        scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0), position = "right") +
        scale_color_manual(values=c(white="white", black="black")) +
        guides(fill=guide_colourbar(title=NULL, barheight=unit(0.65,"npc"), label.position = "left"), color="none") +
        xlab(NULL) + ylab(NULL)
    # if (!missing(title)) {
    #     ggp <- ggp + ggtitle(title)
    # }

    return (ggp)
}

pdf("./overview/mantel_tests_BS.pdf", 5.2, 2.1)
print(manteltest_plot(mt_intra$C, t(mt_intra$P), mt_intra$Cil, mt_intra$Ciu, title="Intra-individual"))
#print(manteltest_plot(mt_inter_aod$C, t(mt_inter_aod$P), mt_inter_aod$Cil, mt_inter_aod$Ciu, title="Inter-individual (AOD)"))
print(manteltest_plot(mt_inter_doa$C, t(mt_inter_doa$P), mt_inter_doa$Cil, mt_inter_doa$Ciu, title="Inter-individual (DOA)"))
dev.off()

pdf("./overview/mantel_tests_fig1.pdf", 2.9 + 0.122047-.35, 1.4*.92)
# print(manteltest_plot(pmax(mt_intra$C, t(mt_inter_aod$C), na.rm=T),
#                       pmax(mt_intra$P, t(mt_inter_aod$P), na.rm=T),
#                       stars=T, title="AoD"))
print(manteltest_plot(pmax(mt_intra$C, t(mt_inter_doa$C), na.rm=T),
                      pmax(mt_intra$P, t(mt_inter_doa$P), na.rm=T),
                      stars=T))
dev.off()
