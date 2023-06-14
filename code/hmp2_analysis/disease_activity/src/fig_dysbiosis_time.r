
library(dplyr)
library(reshape2)
library(ggplot2)

#source("./common/merge_metadata.r")
source("./common/load_bugs.r")
source("./common/disease_activity.r")
source("./common/theme_nature.r")
source("./common/bug_colors.r")

nactive <- sapply(split(bugs.pcl$meta$active, bugs.pcl$meta$subject), sum)
nsamples <- sapply(split(bugs.pcl$meta$active, bugs.pcl$meta$subject), length)
I <- order(nactive + nsamples * 0.001, decreasing = T)
nactive <- nactive[I]
nsamples <- nsamples[I]

ggplot(data=bugs.pcl$meta, aes(x=week_num, y=factor(subject, levels=names(nactive)))) +
    geom_point(aes(color=active)) +
    facet_grid(diagnosis ~ .) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"))

pdf("./overview/active_disease_counts_across_diagnoses.pdf", 2.433, 1.141)
ggplot(data=rbind(
        cbind(data.frame(n=nactive, typ="Dysbiotic"), bugs.pcl$meta[match(names(nactive), bugs.pcl$meta$subject),]),
        cbind(data.frame(n=nsamples-nactive, typ="Non-dysbiotic"), bugs.pcl$meta[match(names(nactive), bugs.pcl$meta$subject),]))) +
    geom_bar(aes(x=factor(subject, levels=names(nactive)), y=n, fill=factor(typ, levels=c("Non-dysbiotic", "Dysbiotic"))), stat="identity", color="black", size=0.05) +
    facet_grid(. ~ diagnosis, scales = "free_x", space="free_x") +
    scale_fill_manual(values=c(Dysbiotic="red", "Non-dysbiotic"="dodgerblue"), name=NULL) +
    ylab("Count") + xlab("Participant") + scale_y_continuous(expand=c(0, 0)) +
    theme_nature() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
    guides(fill="none")
dev.off()

time_subject_active_plot <- function(df, title = ''){
  #df <- nonIBD_df
  df %>% group_by(subject) %>% tally(active, sort = T) -> df2
  #merge(df, df2, by= 'subject') -> df3
  #df4 <- df3[order(-df3$n),]
  df$subject <- factor(df$subject, levels=df2$subject)
  the_plot <- ggplot(df, aes(week_num, subject, fill = active)) +
    geom_point( alpha = .9, shape = 21, size = 1, stroke = 0.01) +
    scale_fill_manual(breaks = c(TRUE, FALSE),
                       values=c('white', 'maroon'))+
    xlab('Time (week)') +  ylab('Subject') + labs(title = title )+
    theme_nature()+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(the_plot)
}
df_all <- bugs.pcl$meta[, c('subject', 'diagnosis', 'week_num', 'activity_index', 'active' )]
nonIBD_df <- df_all[df_all$diagnosis == 'nonIBD',]
UC_df <- df_all[df_all$diagnosis == 'UC',]
CD_df <- df_all[df_all$diagnosis == 'CD',]
nonIBD_plot <- time_subject_active_plot(nonIBD_df, 'nonIBD') + theme(legend.position="none")
UC_plot <- time_subject_active_plot(UC_df, 'UC') + theme(axis.title.y=element_blank(), legend.position="none")
CD_plot <- time_subject_active_plot(CD_df, 'CD') + theme(axis.title.y=element_blank(), legend.position="none")



df_all %>% group_by(subject) %>% tally(active, sort = T) -> df2
merge(df_all, df2, by= 'subject') -> df3
df4 <- df3[order(-df3$n),]
distribution_plot <- ggplot(df4, aes(x=n)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.5, fill="maroon") +
  xlab('Number of dybiosis') +  ylab('Density') + labs(title = '' )+
  theme_nature()
active_status_plot <- plot_grid(nonIBD_plot,
                                UC_plot,
                                CD_plot,
                                distribution_plot,
                           labels=c("a", "b", "c", "d"),
                           ncol = 4, nrow = 1, align = 'v',
                           hjust=-0.1)
ggsave(filename='./disease_activity/active_status_plot.pdf',
       plot=active_status_plot, width = 183, height = 60, units = "mm", dpi = 350)
