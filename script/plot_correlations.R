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
require(cowplot)

plot_corr <- function(df, graphname, xlab, ylab){
   print (paste(graphname, "correlation:", cor(df$param1, df$param2)))
   textsize <- 7
   p <-  ggplot() +
     geom_point(data=df, aes(x=param1, y=param2, size=param_freq), color='black', alpha=0.4, pch=16) +
     scale_size_continuous(limits = c(-6, -1), range = c(0.1, 1), breaks = c(-5, -4, -3, -2, -1),
                           labels = c(expression(bold('10'^'-5')), expression(bold('10'^'-4')),
                                      expression(bold('10'^'-3')), expression(bold('10'^'-2')),
                                      expression(bold('10'^'-1')))) +
     theme_cowplot(12) +
     theme(plot.title=element_text(size=textsize+1,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
           axis.text=element_text(size=textsize,face="bold",colour = 'black'),
           axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
           axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
           axis.title.x=element_text(size=textsize,face="bold"),
           axis.title.y=element_text(size=textsize,face="bold"),
           legend.position = "right",
           legend.title    = element_text(size=textsize,face="bold"),
           legend.text=element_text(size=textsize,face="bold"),
           legend.justification='center',
           legend.key.size = unit(0.3,"line")) +
     ylab(ylab) +
     xlab(xlab) +
     labs(size="Occurrence\nfreq")
   ggsave(graphname,p,width=2, height=1.5, bg='white', dpi=600)
 }

plot_rep_corr <- function(df, ab_name, graphname){
  print (paste(graphname, "correlation:", cor(df$param1, df$param2)))
  textsize <- 7
  p <-  ggplot() +
    geom_point(data=df, aes(x=param1, y=param2, size=param_freq), color='black', alpha=0.4, pch=16) +
    scale_size_continuous(limits = c(-6, -1), range = c(0.1, 1), breaks = c(-5, -4, -3, -2, -1),
                          labels = c(expression(bold('10'^'-5')), expression(bold('10'^'-4')),
                                     expression(bold('10'^'-3')), expression(bold('10'^'-2')),
                                     expression(bold('10'^'-1')))) +
    ggtitle(ab_name) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize+1,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
          axis.text=element_text(size=textsize,face="bold",colour = 'black'),
          axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
          axis.title.x=element_text(size=textsize,face="bold"),
          axis.title.y=element_text(size=textsize,face="bold"),
          legend.position = "right",
          legend.title    = element_text(size=textsize,face="bold"),
          legend.text=element_text(size=textsize,face="bold"),
          legend.justification='center',
          legend.key.size = unit(0.3,"line")) +
    ylab("Replicate 2") +
    xlab("Replicate 1") +
    labs(size="Occurrence\nfreq")
  ggsave(graphname,p,width=2, height=1.5, bg='white',dpi=600)
}

##### MAIN #####

df <- read_tsv('result/S2HR1_scores_common.tsv') %>%
        filter(avg_freq > 0.00002)

print (length(df$avg_freq))

df_plot <- df %>%
             mutate(param1 = A107_score_rep1) %>%
             mutate(param2 = A107_score_rep2) %>%
             mutate(param_freq = log10(avg_freq))
plot_rep_corr(df_plot, 'COVA1-07', 'graph/cor_rep_COVA107.png')

df_plot <- df %>%
             mutate(param1 = A214_score_rep1) %>%
             mutate(param2 = A214_score_rep2) %>%
             mutate(param_freq = log10(avg_freq))
plot_rep_corr(df_plot, 'COVA2-14', 'graph/cor_rep_COVA214.png')

df_plot <- df %>%
             mutate(param1 = A218_score_rep1) %>%
             mutate(param2 = A218_score_rep2) %>%
             mutate(param_freq = log10(avg_freq))
plot_rep_corr(df_plot, 'COVA2-18', 'graph/cor_rep_COVA218.png')

df_plot <- df %>%
             mutate(param1 = A107_score) %>%
             mutate(param2 = A214_score) %>%
             mutate(param_freq = log10(avg_freq))
plot_corr(df_plot, 'graph/cor_bind_COVA107_vs_COVA214.png', 'COVA1-07 binding score', 'COVA2-14 binding score')

df_plot <- df %>%
             mutate(param1 = A107_score) %>%
             mutate(param2 = A218_score) %>%
             mutate(param_freq = log10(avg_freq))
plot_corr(df_plot, 'graph/cor_bind_COVA107_vs_COVA218.png', 'COVA1-07 binding score', 'COVA2-18 binding score')

df_plot <- df %>%
             mutate(param1 = A214_score) %>%
             mutate(param2 = A218_score) %>%
             mutate(param_freq = log10(avg_freq))
plot_corr(df_plot, 'graph/cor_bind_COVA214_vs_COVA218.png', 'COVA2-14 binding score', 'COVA2-18 binding score')
