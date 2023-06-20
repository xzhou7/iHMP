library(patchwork)
library(ggstatsplot)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP2_HMP2_Core/")
source("./Analysis/Analysis Script/sxt.tools.R")

#prevalence core figure
table(colSums(ST.Pr) == 0)
table(colSums(SK.Pr) == 0)
table(colSums(OR.Pr) == 0)
table(colSums(NS.Pr) == 0)

ST.Pr.Char <- ST.Pr[,colSums(ST.Pr) != 0]
ST.Pr.Char[ST.Pr.Char >= 0.8] <- "core"
ST.Pr.Char[ST.Pr.Char < 0.8 & ST.Pr.Char> 0.2] <- "interm"
ST.Pr.Char[ ST.Pr.Char<= 0.2] <- "oppor"

dat.ST <- data.frame()
dat.ST[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST)[1] <- "Var1"
dim(ST.Pr.Char)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char)[i]
  dat.ST <- merge(dat.ST,dat1, by="Var1", all=T)
}

dat.ST[is.na(dat.ST)] <- 0

dat.ST.melt <- melt(dat.ST,id.vars = "Var1")
listst <- as.character(dat.ST.melt$variable[dat.ST.melt$value==55&dat.ST.melt$Var1=="oppor"])
dat.ST.melt$cater <- "cat1"
dat.ST.melt$cater[which(dat.ST.melt$variable %in% listst)] <- "cat2"

order.ST <- data.frame(t(dat.ST))[2:361,]
order.ST$X1 <- as.numeric(order.ST$X1)
order.ST$X2 <- as.numeric(order.ST$X2)
order.ST$X3 <- as.numeric(order.ST$X3)
order.ST <- order.ST[with(order.ST, order(X1, X2, X3)),]
order.ST$genus <- row.names(order.ST)
order.ST$order <- 1:360

dat.ST.melt$order <- order.ST$order[match(dat.ST.melt$variable,order.ST$genus)]

pcore.percent1 <- ggplot(dat.ST.melt, mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent1 <- pcore.percent1 + facet_wrap(.~cater, scales = "free_y")
pcore.percent1

#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.ST.pdf",pcore.percent1, width=12, height=24, dpi=300)

#seperate IRIS
ST.Pr.Char.IS <- ST.Pr.Char[row.names(ST.Pr.Char) %in% subset(sc,IRIS=="IS")$SubjectID,]
ST.Pr.Char.IR <- ST.Pr.Char[row.names(ST.Pr.Char) %in% subset(sc,IRIS=="IR")$SubjectID,]

dat.ST.IS <- data.frame()
dat.ST.IS[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST.IS)[1] <- "Var1"
dim(ST.Pr.Char.IS)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char.IS[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char.IS)[i]
  dat.ST.IS <- merge(dat.ST.IS,dat1, by="Var1", all=T)
}
dat.ST.IS[is.na(dat.ST.IS)] <- 0
dat.ST.IS$IRIS <- "IS"

dat.ST.IR <- data.frame()
dat.ST.IR[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST.IR)[1] <- "Var1"
dim(ST.Pr.Char.IR)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char.IR[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char.IR)[i]
  dat.ST.IR <- merge(dat.ST.IR,dat1, by="Var1", all=T)
}
dat.ST.IR[is.na(dat.ST.IR)] <- 0
dat.ST.IR$IRIS <- "IR"

#combine two dataset
dat.ST.melt.IRIS <- melt(rbind(dat.ST.IS,dat.ST.IR),id.vars = c("Var1", "IRIS"))
dat.ST.melt.IRIS$cater <- "cat1"
#put all opportunistic microbes as a seperate category
dat.ST.melt.IRIS$cater[which(dat.ST.melt.IRIS$variable %in% listst)] <- "cat2"

#adjust the plot order
orderIS <- data.frame(t(dat.ST.IS))[2:361,]
orderIS$X1 <- as.numeric(orderIS$X1)
orderIS$X2 <- as.numeric(orderIS$X2)
orderIS$X3 <- as.numeric(orderIS$X3)
orderIS <- orderIS[with(orderIS, order(X1, X2, X3)),]
orderIS$genus <- row.names(orderIS)
orderIS$order <- 1:360

#dat.ST.melt.IRIS$order <- dat.ST.melt.IRIS$variable

dat.ST.melt.IRIS$order <- orderIS$order[match(dat.ST.melt.IRIS$variable,orderIS$genus)]
dat.ST.melt.IRIS$IRIS <- factor(dat.ST.melt.IRIS$IRIS, levels=c("IS", "IR"))

pcore.percent3 <- filter(dat.ST.melt.IRIS, cater=="cat1") %>% ggplot(mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent3 <- pcore.percent3 + facet_wrap(.~IRIS)
pcore.percent3

#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.ST.IRIS.pdf",pcore.percent3, width=12, height=24, dpi=300)

subset(dat.ST.melt.IRIS, variable == "Coprococcus")$value[1:3]
chisq.st.result <- data.frame()
for (i in 1:360){
  print(i)
  taxa.i <- as.character(unique(dat.ST.melt.IRIS$variable))[i]
  print(taxa.i)
  
  taxa.st.chi <- chisq.test(cbind(c(subset(dat.ST.melt.IRIS, variable == taxa.i)$value[1],sum(subset(dat.ST.melt.IRIS, variable == taxa.i)$value[2:3])),
                                  c(subset(dat.ST.melt.IRIS, variable == taxa.i)$value[4],sum(subset(dat.ST.melt.IRIS, variable == taxa.i)$value[5:6]))))
  chisq.st.result[i,1] <- taxa.i
  chisq.st.result[i,2] <- taxa.st.chi$p.value
  chisq.st.result[i,3] <- taxa.st.chi$statistic
}

chisq.st.result <- chisq.st.result[!is.na(chisq.st.result$V2),]
colnames(chisq.st.result) <- c("taxa", "p", "X2")
filter(chisq.st.result, p < 0.05)

write.csv(file = "./Analysis/tables/chise.st.csv",chisq.st.result)

subset(dat.ST.melt.IRIS, variable == "Butyrivibrio")$value

taxa.chi <- chisq.test(cbind(subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[1:3],
                             subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[4:6]))
taxa.chi

taxa.fisher <- fisher.test(cbind(subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[1:3],
                                 subset(dat.ST.melt.IRIS, variable == "Barnesiella")$value[4:6]))

#Coprococcus
fisher.test(cbind(c(20,1),c(13,7)))
chisq.test(cbind(c(20,1),c(13,7)))

#Parasutterella
fisher.test(cbind(c(14,7),c(6,14)))
chisq.test(cbind(c(14,7),c(6,14)))

#Intestinimonas
fisher.test(cbind(c(10,11),c(8,12)))
chisq.test(cbind(c(10,11),c(8,12)))

#Butyricimonas
fisher.test(cbind(c(10,11),c(5,15)))
chisq.test(cbind(c(10,11),c(5,15)))

#Butyricicoccus
fisher.test(cbind(c(9,12),c(2,18)))
chisq.test(cbind(c(9,12),c(2,18)))

#Butyrivibrio
fisher.test(cbind(c(4,17),c(0,20)))
chisq.test(cbind(c(4,17),c(0,20)))



##############################################################################################################################################################
ST.Pr.Char.OB <- ST.Pr.Char[row.names(ST.Pr.Char) %in% subset(sc,BMI > 30)$SubjectID,]
ST.Pr.Char.NC <- ST.Pr.Char[row.names(ST.Pr.Char) %in% subset(sc,BMI < 30)$SubjectID,]

dat.ST.OB <- data.frame()
dat.ST.OB[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST.OB)[1] <- "Var1"
dim(ST.Pr.Char.OB)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char.OB[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char.OB)[i]
  dat.ST.OB <- merge(dat.ST.OB,dat1, by="Var1", all=T)
}
dat.ST.OB[is.na(dat.ST.OB)] <- 0
dat.ST.OB$OBSE <- "OB"

dat.ST.NC <- data.frame()
dat.ST.NC[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.ST.NC)[1] <- "Var1"
dim(ST.Pr.Char.NC)
for (i in 1:360){
  print(i)
  dat1 <- data.frame(table(sapply(ST.Pr.Char.NC[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(ST.Pr.Char.NC)[i]
  dat.ST.NC <- merge(dat.ST.NC,dat1, by="Var1", all=T)
}
dat.ST.NC[is.na(dat.ST.NC)] <- 0
dat.ST.NC$OBSE <- "NC"

#combine two dataset
dat.ST.melt.OBSE <- melt(rbind(dat.ST.OB,dat.ST.NC),id.vars = c("Var1", "OBSE"))
dat.ST.melt.OBSE$cater <- "cat1"
#put all opportunistic microbes as a seperate category
dat.ST.melt.OBSE$cater[which(dat.ST.melt.OBSE$variable %in% listst)] <- "cat2"

#adjust the plot order
orderNC <- data.frame(t(dat.ST.NC))[2:361,]
orderNC$X1 <- as.numeric(orderNC$X1)
orderNC$X2 <- as.numeric(orderNC$X2)
orderNC$X3 <- as.numeric(orderNC$X3)
orderNC <- orderNC[with(orderNC, order(X1, X2, X3)),]
orderNC$genus <- row.names(orderNC)
orderNC$order <- 1:360

#dat.ST.melt.IRIS$order <- dat.ST.melt.IRIS$variable

dat.ST.melt.OBSE$order <- orderNC$order[match(dat.ST.melt.OBSE$variable,orderNC$genus)]
dat.ST.melt.OBSE$OBSE <- factor(dat.ST.melt.OBSE$OBSE, levels=c("NC", "OB"))

pcore.percent3.1 <- filter(dat.ST.melt.OBSE, cater=="cat1") %>% ggplot(mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent3.1 <- pcore.percent3.1 + facet_grid(.~OBSE, scales = "free")
pcore.percent3.1
#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.ST.OBSE.pdf",pcore.percent3.1, width=12, height=24, dpi=300)

chisq.st.ob.result <- data.frame()
for (i in 1:360){
  print(i)
  taxa.i <- as.character(unique(dat.ST.melt.OBSE$variable))[i]
  print(taxa.i)
  
  taxa.st.chi <- chisq.test(cbind(c(subset(dat.ST.melt.OBSE, variable == taxa.i)$value[1],sum(subset(dat.ST.melt.OBSE, variable == taxa.i)$value[2:3])),
                                  c(subset(dat.ST.melt.OBSE, variable == taxa.i)$value[4],sum(subset(dat.ST.melt.OBSE, variable == taxa.i)$value[5:6]))))
  chisq.st.ob.result[i,1] <- taxa.i
  chisq.st.ob.result[i,2] <- taxa.st.chi$p.value
  chisq.st.ob.result[i,3] <- taxa.st.chi$statistic
}

chisq.st.ob.result <- chisq.st.ob.result[!is.na(chisq.st.ob.result$V2),]
colnames(chisq.st.ob.result) <- c("taxa", "p", "X2")
filter(chisq.st.ob.result, p < 0.05)

#write.csv(file = "./Analysis/tables/chise.st.ob.csv",chisq.st.ob.result)
###############################################################################################################################################################
OR.Pr.Char <- OR.Pr[,colSums(OR.Pr) != 0]
OR.Pr.Char[OR.Pr.Char >= 0.8] <- "core"
OR.Pr.Char[OR.Pr.Char < 0.8 & OR.Pr.Char> 0.2] <- "interm"
OR.Pr.Char[ OR.Pr.Char<= 0.2] <- "oppor"

dat.OR <- data.frame()
dat.OR[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.OR)[1] <- "Var1"
dim(OR.Pr.Char)
for (i in 1:392){
  print(i)
  dat1 <- data.frame(table(sapply(OR.Pr.Char[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(OR.Pr.Char)[i]
  dat.OR <- merge(dat.OR,dat1, by="Var1", all=T)
}

dat.OR[is.na(dat.OR)] <- 0

dat.OR.melt <- melt(dat.OR,id.vars = "Var1")
listor <- as.character(dat.OR.melt$variable[dat.OR.melt$value==67&dat.OR.melt$Var1=="oppor"])

OR.Pr.Char.OB <- OR.Pr.Char[row.names(OR.Pr.Char) %in% subset(sc,BMI > 30)$SubjectID,]
OR.Pr.Char.NC <- OR.Pr.Char[row.names(OR.Pr.Char) %in% subset(sc,BMI < 30)$SubjectID,]

dat.OR.OB <- data.frame()
dat.OR.OB[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.OR.OB)[1] <- "Var1"
dim(OR.Pr.Char.OB)
for (i in 1:392){
  print(i)
  dat1 <- data.frame(table(sapply(OR.Pr.Char.OB[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(OR.Pr.Char.OB)[i]
  dat.OR.OB <- merge(dat.OR.OB,dat1, by="Var1", all=T)
}
dat.OR.OB[is.na(dat.OR.OB)] <- 0
dat.OR.OB$OBSE <- "OB"

dat.OR.NC <- data.frame()
dat.OR.NC[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.OR.NC)[1] <- "Var1"
dim(OR.Pr.Char.NC)
for (i in 1:392){
  print(i)
  dat1 <- data.frame(table(sapply(OR.Pr.Char.NC[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(OR.Pr.Char.NC)[i]
  dat.OR.NC <- merge(dat.OR.NC,dat1, by="Var1", all=T)
}
dat.OR.NC[is.na(dat.OR.NC)] <- 0
dat.OR.NC$OBSE <- "NC"

#combine two dataset
dat.OR.melt.OBSE <- melt(rbind(dat.OR.OB,dat.OR.NC),id.vars = c("Var1", "OBSE"))
dat.OR.melt.OBSE$cater <- "cat1"
#put all opportunistic microbes as a seperate category
dat.OR.melt.OBSE$cater[which(dat.OR.melt.OBSE$variable %in% listor)] <- "cat2"

#adjust the plot order
orderNC <- data.frame(t(dat.OR.NC))[2:393,]
orderNC$X1 <- as.numeric(orderNC$X1)
orderNC$X2 <- as.numeric(orderNC$X2)
orderNC$X3 <- as.numeric(orderNC$X3)
orderNC <- orderNC[with(orderNC, order(X1, X2, X3)),]
orderNC$genus <- row.names(orderNC)
orderNC$order <- 1:392

#dat.ST.melt.IRIS$order <- dat.ST.melt.IRIS$variable

dat.OR.melt.OBSE$order <- orderNC$order[match(dat.OR.melt.OBSE$variable,orderNC$genus)]
dat.OR.melt.OBSE$OBSE <- factor(dat.OR.melt.OBSE$OBSE, levels=c("NC", "OB"))

pcore.percent3.2 <- filter(dat.OR.melt.OBSE, cater=="cat1") %>% ggplot(mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent3.2 <- pcore.percent3.2 + facet_grid(.~OBSE, scales = "free")
pcore.percent3.2
#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.OR.OBSE.pdf",pcore.percent3.2, width=12, height=28, dpi=300)

chisq.or.ob.result <- data.frame()
for (i in 1:392){
  print(i)
  taxa.i <- as.character(unique(dat.OR.melt.OBSE$variable))[i]
  print(taxa.i)
  
  taxa.or.chi <- chisq.test(cbind(c(subset(dat.OR.melt.OBSE, variable == taxa.i)$value[1],sum(subset(dat.OR.melt.OBSE, variable == taxa.i)$value[2:3])),
                                  c(subset(dat.OR.melt.OBSE, variable == taxa.i)$value[4],sum(subset(dat.OR.melt.OBSE, variable == taxa.i)$value[5:6]))))
  
  chisq.or.ob.result[i,1] <- taxa.i
  chisq.or.ob.result[i,2] <- taxa.or.chi$p.value
  chisq.or.ob.result[i,3] <- taxa.or.chi$statistic
}

chisq.or.ob.result <- chisq.or.ob.result[!is.na(chisq.or.ob.result$V2),]
colnames(chisq.or.ob.result) <- c("taxa", "p", "X2")
filter(chisq.or.ob.result, p < 0.05)
#write.csv(file = "./Analysis/tables/chise.or.ob.csv",chisq.or.ob.result)

###############################################################################################################################################################
SK.Pr.Char.IS <- SK.Pr.Char[row.names(SK.Pr.Char) %in% subset(sc,IRIS=="IS")$SubjectID,]
SK.Pr.Char.IR <- SK.Pr.Char[row.names(SK.Pr.Char) %in% subset(sc,IRIS=="IR")$SubjectID,]

dat.SK.IS <- data.frame()
dat.SK.IS[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.SK.IS)[1] <- "Var1"
dim(SK.Pr.Char.IS)
for (i in 1:1110){
  print(i)
  dat1 <- data.frame(table(sapply(SK.Pr.Char.IS[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(SK.Pr.Char.IS)[i]
  dat.SK.IS <- merge(dat.SK.IS,dat1, by="Var1", all=T)
}
dat.SK.IS[is.na(dat.SK.IS)] <- 0
dat.SK.IS$IRIS <- "IS"

dat.SK.IR <- data.frame()
dat.SK.IR[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.SK.IR)[1] <- "Var1"
dim(SK.Pr.Char.IR)
for (i in 1:1110){
  print(i)
  dat1 <- data.frame(table(sapply(SK.Pr.Char.IR[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(SK.Pr.Char.IR)[i]
  dat.SK.IR <- merge(dat.SK.IR,dat1, by="Var1", all=T)
}
dat.SK.IR[is.na(dat.SK.IR)] <- 0
dat.SK.IR$IRIS <- "IR"

#combine two dataset
dat.SK.melt.IRIS <- melt(rbind(dat.SK.IS,dat.SK.IR),id.vars = c("Var1", "IRIS"))
dat.SK.melt.IRIS$cater <- "cat1"

dat.SK <- data.frame()
dat.SK[1:3,1] <- c("core", "interm", "oppor")
colnames(dat.SK)[1] <- "Var1"
dim(SK.Pr.Char)
for (i in 1:1110){
  print(i)
  dat1 <- data.frame(table(sapply(SK.Pr.Char[,i], function(x) x)))
  colnames(dat1)[2] <- colnames(SK.Pr.Char)[i]
  dat.SK <- merge(dat.SK,dat1, by="Var1", all=T)
}

dat.SK[is.na(dat.SK)] <- 0

dat.SK.melt <- melt(dat.SK,id.vars = "Var1")
listsk <- as.character(dat.SK.melt$variable[dat.SK.melt$value==69&dat.SK.melt$Var1=="oppor"])

dat.SK.melt.IRIS$cater[which(dat.SK.melt.IRIS$variable %in% listsk)] <- "cat2"

#adjust the plot order
orderIS.sk <- data.frame(t(dat.SK.IS))[2:1111,]
orderIS.sk$X1 <- as.numeric(orderIS.sk$X1)
orderIS.sk$X2 <- as.numeric(orderIS.sk$X2)
orderIS.sk$X3 <- as.numeric(orderIS.sk$X3)
orderIS.sk <- orderIS.sk[with(orderIS.sk, order(X1, X2, X3)),]
orderIS.sk$genus <- row.names(orderIS.sk)
orderIS.sk$order <- 1:1110

#dat.SK.melt.IRIS$order <- dat.SK.melt.IRIS$variable

dat.SK.melt.IRIS$order <- orderIS.sk$order[match(dat.SK.melt.IRIS$variable,orderIS.sk$genus)]
dat.SK.melt.IRIS$IRIS <- factor(dat.SK.melt.IRIS$IRIS, levels=c("IS", "IR"))

pcore.percent3.sk <- filter(dat.SK.melt.IRIS, cater=="cat1") %>% ggplot(mapping=aes(y=reorder(variable,order), x=value, fill=Var1)) + geom_bar(position = "stack",stat="identity")
pcore.percent3.sk <- pcore.percent3.sk + facet_wrap(.~IRIS)
pcore.percent3.sk

#ggsave(filename = "./Analysis/Suppl.figure/Core.microbiome.Pre.SK.IRIS.pdf",pcore.percent3.sk, width=12, height=48, dpi=300)

#chi sq test for all taxa in skin
chisq.sk.result <- data.frame()
for (i in 1:1110){
  print(i)
  taxa.i <- as.character(unique(dat.SK.melt.IRIS$variable))[i]
  print(taxa.i)
  
  taxa.sk.chi <- chisq.test(cbind(c(subset(dat.SK.melt.IRIS, variable == taxa.i)$value[1],sum(subset(dat.SK.melt.IRIS, variable == taxa.i)$value[2:3])),
                                  c(subset(dat.SK.melt.IRIS, variable == taxa.i)$value[4],sum(subset(dat.SK.melt.IRIS, variable == taxa.i)$value[5:6]))))
  chisq.sk.result[i,1] <- taxa.i
  chisq.sk.result[i,2] <- taxa.sk.chi$p.value
  chisq.sk.result[i,3] <- taxa.sk.chi$statistic
}
chisq.sk.result <- chisq.sk.result[!is.na(chisq.sk.result$V2),]
colnames(chisq.sk.result) <- c("taxa", "p", "X2")
filter(chisq.sk.result, p < 0.05)

#write.csv(file = "./Analysis/tables/chise.sk.csv",chisq.sk.result)

#Rothia
cbind(subset(dat.SK.melt.IRIS, variable == "Rothia")$value[1:3],subset(dat.SK.melt.IRIS, variable == "Rothia")$value[4:6])
chisq.test(cbind(c(2,20),c(6,20)))

###############################################################################################################################################################
###############################################################################################################################################################
#core microbiome number by IS and IR individual
stool.cr.IRIS <- ggplot(filter(core.bysubject.meta, IRIS!="Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")),y=observed, fill=IRIS)) + geom_jitter() + geom_boxplot(alpha=0.4, outlier.alpha = 0)
stool.cr.IRIS <- stool.cr.IRIS + base_theme + scale_fill_manual(values = iris_color) + stat_compare_means(method = "t.test")
stool.cr.IRIS <- stool.cr.IRIS + xlab("") + ylab("Number of Core Microbiome") + ggtitle("Stool")
stool.cr.IRIS

skin.cr.IRIS <- ggplot(filter(core.bysubject.meta.sk, IRIS!="Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")),y=observed, fill=IRIS)) + geom_jitter() + geom_boxplot(alpha=0.4, outlier.alpha = 0)
skin.cr.IRIS <- skin.cr.IRIS + base_theme + scale_fill_manual(values = iris_color) + stat_compare_means(method = "t.test")
skin.cr.IRIS <- skin.cr.IRIS + xlab("") + ylab("Number of Core Microbiome") + ggtitle("Skin")
skin.cr.IRIS

oral.cr.IRIS <- ggplot(filter(core.bysubject.meta.or, IRIS!="Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")),y=observed, fill=IRIS)) + geom_jitter() + geom_boxplot(alpha=0.4, outlier.alpha = 0)
oral.cr.IRIS <- oral.cr.IRIS + base_theme + scale_fill_manual(values = iris_color) + stat_compare_means(method = "t.test")
oral.cr.IRIS <- oral.cr.IRIS + xlab("") + ylab("Number of Core Microbiome") + ggtitle("Oral")
oral.cr.IRIS

nasal.cr.IRIS <- ggplot(filter(core.bysubject.meta.ns, IRIS!="Unknown"), aes(x=factor(IRIS, levels=c("IS", "IR")),y=observed, fill=IRIS)) + geom_jitter() + geom_boxplot(alpha=0.4, outlier.alpha = 0)
nasal.cr.IRIS <- nasal.cr.IRIS + base_theme + scale_fill_manual(values = iris_color) + stat_compare_means(method = "t.test")
nasal.cr.IRIS <- nasal.cr.IRIS + xlab("") + ylab("Number of Core Microbiome") + ggtitle("Nasal")
nasal.cr.IRIS

p.cr.IRIS <- stool.cr.IRIS + skin.cr.IRIS + oral.cr.IRIS + nasal.cr.IRIS +  plot_layout(guides = "collect")
p.cr.IRIS

#get t statistics of the comparasion
t.test(filter(core.bysubject.meta, IRIS=="IS")$observed, filter(core.bysubject.meta, IRIS=="IR")$observed)
t.test(filter(core.bysubject.meta.sk, IRIS=="IS")$observed, filter(core.bysubject.meta.sk, IRIS=="IR")$observed)

###############################
load("./Analysis/Robject/DetailedPhyloseq.RData")
st.ns.tax <- tax_table(physeqGenus_ST) %>% data.frame()
sk.or.tax <- tax_table(physeqGenus_SK) %>% data.frame()

dim(st.ns.tax)
dim(sk.or.tax)

#stool
stool.core.taxa <- ST.Pr %>% colMeans() %>% data.frame() %>% merge(st.ns.tax,by.x="row.names", by.y="Genus", all.x=T)
stool.core.taxa$Cate <- "Medi"
stool.core.taxa$Cate[stool.core.taxa$. < 0.2] <- "Oppr"
stool.core.taxa$Cate[stool.core.taxa$. > 0.8] <- "Core"

stool.core.taxa$phylum_core <- "Other"
stool.core.taxa$phylum_core[stool.core.taxa$Phylum == "Actinobacteria"] <- "Actinobacteria"
stool.core.taxa$phylum_core[stool.core.taxa$Phylum == "Bacteroidetes"] <- "Bacteroidetes"
stool.core.taxa$phylum_core[stool.core.taxa$Phylum == "Firmicutes"] <- "Firmicutes"
stool.core.taxa$phylum_core[stool.core.taxa$Phylum == "Proteobacteria"] <- "Proteobacteria"

#skin
skin.core.taxa <- SK.Pr %>% colMeans() %>% data.frame() %>% merge(sk.or.tax,by.x="row.names", by.y="Genus", all.x=T)
skin.core.taxa$Cate <- "Medi"
skin.core.taxa$Cate[skin.core.taxa$. < 0.2] <- "Oppr"
skin.core.taxa$Cate[skin.core.taxa$. > 0.8] <- "Core"

skin.core.taxa$phylum_core <- "Other"
skin.core.taxa$phylum_core[skin.core.taxa$Phylum == "Actinobacteria"] <- "Actinobacteria"
skin.core.taxa$phylum_core[skin.core.taxa$Phylum == "Bacteroidetes"] <- "Bacteroidetes"
skin.core.taxa$phylum_core[skin.core.taxa$Phylum == "Firmicutes"] <- "Firmicutes"
skin.core.taxa$phylum_core[skin.core.taxa$Phylum == "Proteobacteria"] <- "Proteobacteria"

#oral
oral.core.taxa <- OR.Pr %>% colMeans() %>% data.frame() %>% merge(sk.or.tax,by.x="row.names", by.y="Genus", all.x=T)
oral.core.taxa$Cate <- "Medi"
oral.core.taxa$Cate[oral.core.taxa$. < 0.2] <- "Oppr"
oral.core.taxa$Cate[oral.core.taxa$. > 0.8] <- "Core"

oral.core.taxa$phylum_core <- "Other"
oral.core.taxa$phylum_core[oral.core.taxa$Phylum == "Actinobacteria"] <- "Actinobacteria"
oral.core.taxa$phylum_core[oral.core.taxa$Phylum == "Bacteroidetes"] <- "Bacteroidetes"
oral.core.taxa$phylum_core[oral.core.taxa$Phylum == "Firmicutes"] <- "Firmicutes"
oral.core.taxa$phylum_core[oral.core.taxa$Phylum == "Proteobacteria"] <- "Proteobacteria"

#nasal
nasal.core.taxa <- NS.Pr %>% colMeans() %>% data.frame() %>% merge(st.ns.tax,by.x="row.names", by.y="Genus", all.x=T)
nasal.core.taxa$Cate <- "Medi"
nasal.core.taxa$Cate[nasal.core.taxa$. < 0.2] <- "Oppr"
nasal.core.taxa$Cate[nasal.core.taxa$. > 0.8] <- "Core"

nasal.core.taxa$phylum_core <- "Other"
nasal.core.taxa$phylum_core[nasal.core.taxa$Phylum == "Actinobacteria"] <- "Actinobacteria"
nasal.core.taxa$phylum_core[nasal.core.taxa$Phylum == "Bacteroidetes"] <- "Bacteroidetes"
nasal.core.taxa$phylum_core[nasal.core.taxa$Phylum == "Firmicutes"] <- "Firmicutes"
nasal.core.taxa$phylum_core[nasal.core.taxa$Phylum == "Proteobacteria"] <- "Proteobacteria"

percent.core.stool <- stool.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>%
  ggplot(aes(x = factor(Cate), y = perc*100, fill = factor(phylum_core))) +
  geom_bar(stat="identity", width = 0.7) + base_theme + scale_fill_manual(values = phylum_color) + ylab("Percentage of Phylum (%)") + 
  xlab("") + ggtitle("Stool Microbiome")
 
percent.core.stool

percent.core.skin<- skin.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>%
  ggplot(aes(x = factor(Cate), y = perc*100, fill = factor(phylum_core))) +
  geom_bar(stat="identity", width = 0.7) + base_theme + scale_fill_manual(values = phylum_color) + ylab("Percentage of Phylum (%)") + 
  xlab("") + ggtitle("Skin Microbiome")

percent.core.skin

percent.core.oral<- oral.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>%
  ggplot(aes(x = factor(Cate), y = perc*100, fill = factor(phylum_core))) +
  geom_bar(stat="identity", width = 0.7) + base_theme + scale_fill_manual(values = phylum_color) + ylab("Percentage of Phylum (%)") + 
  xlab("") + ggtitle("Oral Microbiome")

percent.core.oral

percent.core.nasal<- nasal.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>%
  ggplot(aes(x = factor(Cate), y = perc*100, fill = factor(phylum_core))) +
  geom_bar(stat="identity", width = 0.7) + base_theme + scale_fill_manual(values = phylum_color) + ylab("Percentage of Phylum (%)") + 
  xlab("") + ggtitle("Nasal Microbiome")

percent.core.nasal

p.percent.core <- percent.core.stool + percent.core.skin + percent.core.oral + percent.core.nasal + plot_layout(guides = "collect")
p.percent.core
#ggsave(filename = "./Analysis/Suppl.figure/p.percent.core.pdf", p.percent.core, width = 10, height = 7, dpi = 300)


p.stool <- ggbarstats(
  data= stool.core.taxa,
  x= phylum_core,
  y= Cate,
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  title = "Stool Microbiome By Percent"
)
p.stool <- p.stool + scale_fill_manual(values = phylum_color)
p.stool
#ggsave(filename = "./Analysis/Suppl.figure/p.percent.core.stool.2.pdf", p.stool, width = 10, height = 7, dpi = 300)

p.skin <- ggbarstats(
  data= skin.core.taxa,
  x= phylum_core,
  y= Cate,
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  title = "Skin Microbiome By Percent"
)
p.skin <- p.skin + scale_fill_manual(values = phylum_color)
p.skin
#ggsave(filename = "./Analysis/Suppl.figure/p.percent.core.skin.2.pdf", p.skin, width = 10, height = 7, dpi = 300)

p.oral <- ggbarstats(
  data= oral.core.taxa,
  x= phylum_core,
  y= Cate,
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  title = "Oral Microbiome By Percent"
)
p.oral <- p.oral + scale_fill_manual(values = phylum_color)
p.oral
#ggsave(filename = "./Analysis/Suppl.figure/p.percent.core.oral.2.pdf", p.oral, width = 10, height = 7, dpi = 300)

p.nasal <- ggbarstats(
  data= nasal.core.taxa,
  x= phylum_core,
  y= Cate,
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  title = "Nasal Microbiome By Percent"
)
p.nasal <- p.nasal + scale_fill_manual(values = phylum_color)
p.nasal
#ggsave(filename = "./Analysis/Suppl.figure/p.percent.core.nasal.2.pdf", p.nasal, width = 10, height = 7, dpi = 300)

percent.core.stool.table <- stool.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
percent.core.stool.table$bodysite <- "Stool"

percent.core.skin.table <- skin.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
percent.core.skin.table$bodysite <- "Skin"

percent.core.oral.table <- oral.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 
percent.core.oral.table$bodysite <- "Oral"

percent.core.nasal.table <- nasal.core.taxa %>% group_by(Cate, phylum_core) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
percent.core.nasal.table$bodysite <- "Nasal"


percent.core <- rbind(percent.core.stool.table,
                      percent.core.skin.table,
                      percent.core.oral.table,
                      percent.core.nasal.table)
#write.csv(file = "./Analysis/tables/Core_PercentPhylum.csv",percent.core)
percent.core



