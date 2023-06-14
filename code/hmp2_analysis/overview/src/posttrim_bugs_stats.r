
source("./common/load_bugs.r")
source("./common/load_bugs_posttrim.r")

species_pretrim.pcl <- pcl.only(bugs.pcl, rank="s")
species_posttrim.pcl <- pcl.only(bugs_posttrim.pcl, rank="s")

# Create matched tables
cmn_species <- union(colnames(species_posttrim.pcl$x), colnames(species_pretrim.pcl$x))
species_pretrim.pcl <- species_pretrim.pcl %>% pcl.reorder.f(cmn_species, na.val=0)
species_posttrim.pcl <- species_posttrim.pcl %>% pcl.reorder.f(cmn_species, na.val=0)

# Species gain/loss overall stats
gained <- pcl.apply.s(species_posttrim.pcl, sum((x>0) & (species_pretrim.pcl$x[Name,]==0)))
summary(gained)

lost <- pcl.apply.s(species_posttrim.pcl, sum((x==0) & (species_pretrim.pcl$x[Name,]>0)))
summary(lost)

common <- pcl.apply.s(species_posttrim.pcl, sum((x>0) & (species_pretrim.pcl$x[Name,]>0)))
summary(common)

jaccard_sim <- common / (common + gained + lost)
summary(jaccard_sim)

BC <- pcl.apply.s(species_posttrim.pcl, {
    sum(abs(species_pretrim.pcl$x[Name,] - x)) / 2
})
summary(BC)

# Frequently gained/lost species
post_to_pre <- match(rownames(species_posttrim.pcl$x), rownames(species_pretrim.pcl$x))
sp_gained <- pcl.apply.f(species_posttrim.pcl, sum((x>0) & (species_pretrim.pcl$x[post_to_pre,Name]==0)))
sp_lost <- pcl.apply.f(species_posttrim.pcl, sum((x==0) & (species_pretrim.pcl$x[post_to_pre,Name]>0)))

# Greatest gain/loss of abundance
sp_gain_ab <- pcl.apply.f(species_posttrim.pcl, sum(pmax(0, x - species_pretrim.pcl$x[post_to_pre,Name])))
sp_loss_ab <- pcl.apply.f(species_posttrim.pcl, -sum(pmin(0, x - species_pretrim.pcl$x[post_to_pre,Name])))
