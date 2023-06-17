#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggrepel)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_bind_vs_exp <- function(df, df_special, df_others, df_WT, param, title, graphname){
  textsize <- 7
  print (paste(graphname, "correlation: ", cor(df$exp_score, df[[param]])))
  palette <- c('grey',qualpal(n = 5, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  palette <- c(palette[3], palette[2], palette[4], palette[1], palette[5])
  p <-  ggplot() +
          geom_point(data=df_others, aes(x=exp_score, y=.data[[param]], size=log10(avg_freq), color='Others'), alpha=0.4, pch=16) +
          geom_point(data=df_WT, aes(x=exp_score, y=.data[[param]], size=log10(avg_freq), color='WT'), alpha=0.7, pch=16) +
          geom_point(data=df_special, aes(x=exp_score, y=.data[[param]], size=log10(avg_freq), color=grouping), alpha=0.7, pch=16) +
          scale_color_manual('', values=palette,drop=FALSE) +
          geom_text_repel(data=df_special, aes(x=exp_score, y=.data[[param]],label=mut),
                          color="black", min.segment.length=0, segment.size=0.3, size=2, force=50, force_pull=0,
                          seed=6, max.overlaps = Inf) +
          scale_size_continuous(limits = c(-6, -1), range = c(0.1, 1.5), breaks = c(-5, -4, -3, -2, -1),
                                labels = c(expression(bold('10'^'-5')), expression(bold('10'^'-4')),
                                           expression(bold('10'^'-3')), expression(bold('10'^'-2')),
                                           expression(bold('10'^'-1')))) +
          ggtitle(title) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold", margin=margin(t=1,b=1)),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab("Binding score") +
          xlab("Expression score") +
          labs(size="Avg freq")
  ggsave(graphname,p,width=3, height=2, bg='white')
  }

coloring <- function(mut){
  if (mut %in% c("T961F", "V987C", "Q1010W", "Q1005R")){return ("Low binding")}
  else if (mut %in% c("V915H")){return ("High binding")}
  else if (mut %in% c("WT")){return ("WT")}
  else if (mut %in% c("D950N", "Q954H", "N969K", "L981F", "S982A", "T1027I")){return ("Natural variants")}
  else {return ('Others')}
  } 

##### MAIN #####
group_level <- c('High binding', 'Low binding', 'Natural variants')

df <- read_tsv('result/S2HR1_scores_common.tsv') %>%
        filter(avg_freq > 0.00002) %>%
        mutate(grouping=mapply(coloring,mut)) %>%
        filter(mut_class %in% c('missense', 'WT'))
df_special <- filter(df, grouping %in% group_level)
df_WT <- filter(df, grouping=='WT')
df_others <- filter(df, grouping=='Others')
plot_bind_vs_exp(df, df_special, df_others, df_WT, 'A107_score', 'COVA1-07', 'graph/exp_vs_bind_COVA107.png')
plot_bind_vs_exp(df, df_special, df_others, df_WT, 'A214_score', 'COVA2-14', 'graph/exp_vs_bind_COVA214.png')
plot_bind_vs_exp(df, df_special, df_others, df_WT, 'A218_score', 'COVA2-18', 'graph/exp_vs_bind_COVA218.png')
