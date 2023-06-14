
source("./common/load_bugs.r")

library(ggplot2)

bugs.pcl.nn <- pcl.nicenames(bugs.pcl)

df <- cbind(bugs.pcl.nn$meta, as.data.frame(bugs.pcl.nn$x))
df <- df[!is.na(df$consent_age),]

pdf("./overview/all_bugs_vs_age.pdf", 5, 4)
for (c in colnames(df)) {
    if (is.numeric(df[,c]) && mean(is.na(df[,c])) < 0.05 && mean(na.omit(df[,c])>0) > 0.05) {
        print(ggplot(data=df, aes_string(x="consent_age", y=sprintf("`%s`", c))) +
            geom_smooth() + geom_point())
        print(ggplot(data=df, aes_string(x="adult", y=sprintf("`%s`", c))) +
            geom_boxplot())
    }
}
dev.off()

