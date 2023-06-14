
source("./common/disease_colors.r")

library(viridis)
make_heatmap <- function(pcl, ...) {
    pcl.heatmap(pcl,
                meta=c("diagnosis"), annotation_colors=list(diagnosis=hmp2_disease_colors),
                show_rownames=T, color=viridis(100), logspace=T, minx=1e-4, zerospecial=NA, ...)
}

# --------------------------------------------
# Metabolites

source("./common/load_metabolites.r")

metabolites.pcoa.ord <- pcl.pcoa(metabolites.pcl.nrm)
pdf("./overview/metabolites_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(metabolites.pcl, metabolites.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

metabolites.tsne.ord <- pcl.tsne(metabolites.pcl.nrm)
pdf("./overview/metabolites_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(metabolites.pcl, metabolites.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

metabolites.pcl.nn <- metabolites.pcl.nrm
colnames(metabolites.pcl.nn$x) <- metabolites.pcl.nn$metaf$Metabolite

pdf("./overview/metabolites_heatmap_top50_ab.pdf", 18, 12, onefile=F)
make_heatmap(metabolites.pcl.nn %>% pcl.top.f(mean(x), n=50))
dev.off()

pdf("./overview/metabolites_heatmap_top50_ab_nohclust.pdf", 18, 12, onefile=F)
make_heatmap(metabolites.pcl.nn %>% pcl.top.f(mean(x), n=50) %>%
                 pcl.sort.s(sprintf("%g%s", as.numeric(diagnosis), site_sub_coll)), cluster_cols=F)
dev.off()


# --------------------------------------------
# Proteins

source("./common/load_proteins.r")

proteins.pcl.nrm <- pcl.normalize(proteins.pcl)
proteins.pcoa.ord <- pcl.pcoa(proteins.pcl.nrm)
pdf("./overview/proteins_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(proteins.pcl.nrm, proteins.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

proteins.tsne.ord <- pcl.tsne(proteins.pcl.nrm)
pdf("./overview/proteins_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(proteins.pcl.nrm, proteins.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

proteins.pcl.bn <- proteins.pcl
proteins.pcl.bn$x <- proteins.pcl.bn$x > 0
proteins.pcoa.ord <- pcl.pcoa(proteins.pcl.bn)
pdf("./overview/proteins_pcoa_jacc.pdf", 5.5, 5)
pcl.ordplot(proteins.pcl.bn, proteins.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

proteins.tsne.ord <- pcl.tsne(proteins.pcl.bn)
pdf("./overview/proteins_tsne_jacc.pdf", 5.5, 5)
pcl.ordplot(proteins.pcl.bn, proteins.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

pdf("./overview/proteins_heatmap_top50_ab.pdf", 18, 12, onefile=F)
make_heatmap(proteins.pcl %>% pcl.top.f(mean(x), n=50))
dev.off()

pdf("./overview/proteins_heatmap_top50_ab_nohclust.pdf", 18, 12, onefile=F)
make_heatmap(proteins.pcl %>% pcl.top.f(mean(x), n=50) %>%
                 pcl.sort.s(sprintf("%g%s", as.numeric(diagnosis), site_sub_coll)), cluster_cols=F)
dev.off()


# --------------------------------------------
# Viruses

source("./common/load_bugs.r")

viruses.clean.pcl <- viruses.pcl %>%
    pcl.filter.s(any(x>0)) %>%
    pcl.normalize %>%
    pcl.nicenames

viruses.pcoa.ord <- pcl.pcoa(viruses.clean.pcl)
pdf("./overview/viruses_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(viruses.pcl, viruses.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

viruses.tsne.ord <- pcl.tsne(viruses.clean.pcl)
pdf("./overview/viruses_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(viruses.clean.pcl, viruses.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

pdf("./overview/viruses_heatmap_top50_ab.pdf", 18, 12, onefile=F)
make_heatmap(viruses.clean.pcl %>% pcl.top.f(mean(x), n=10))
dev.off()

pdf("./overview/viruses_heatmap_top50_ab_nohclust.pdf", 18, 12, onefile=F)
make_heatmap(viruses.clean.pcl %>% pcl.top.f(mean(x), n=10) %>%
                 pcl.sort.s(sprintf("%g%s", as.numeric(diagnosis), site_sub_coll)), cluster_cols=F)
dev.off()


# --------------------------------------------
# MTX

source("./common/load_kos.r")
source("./common/load_ecs.r")

ko.rna.unstrat.pcoa.ord <- pcl.pcoa(ko.rna.unstrat.pcl %>% pcl.normalize)
pdf("./overview/mtx_ko_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(ko.rna.unstrat.pcl, ko.rna.unstrat.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

ko.bugs.rna.pcoa.ord <- pcl.pcoa(bugs.fromko.rna.pcl %>% pcl.normalize)
pdf("./overview/mtx_kobug_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(bugs.fromko.rna.pcl, ko.bugs.rna.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

ec.rna.unstrat.pcoa.ord <- pcl.pcoa(ec.rna.unstrat.pcl %>% pcl.normalize)
pdf("./overview/mtx_ec_pcoa_bc.pdf", 5.5, 5)
pcl.ordplot(ec.rna.unstrat.pcl, ec.rna.unstrat.pcoa.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

ko.rna.unstrat.tsne.ord <- pcl.tsne(ko.rna.unstrat.pcl %>% pcl.normalize)
pdf("./overview/mtx_ko_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(ko.rna.unstrat.pcl, ko.rna.unstrat.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

ko.bugs.rna.tsne.ord <- pcl.tsne(bugs.fromko.rna.pcl %>% pcl.normalize)
pdf("./overview/mtx_kobug_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(bugs.fromko.rna.pcl, ko.bugs.rna.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

ec.rna.unstrat.tsne.ord <- pcl.tsne(ec.rna.unstrat.pcl %>% pcl.normalize)
pdf("./overview/mtx_ec_tsne_bc.pdf", 5.5, 5)
pcl.ordplot(ec.rna.unstrat.pcl, ec.rna.unstrat.tsne.ord,
            colour="diagnosis", colour_override=hmp2_disease_colors)
dev.off()

