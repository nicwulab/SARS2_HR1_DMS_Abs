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
library(sinaplot)
library(ggforce)
require(cowplot)

t_test <- function(df, class_1, class_2){
  p_value_A107 <- t.test(filter(df, mut_class==class_1)$A107_score, filter(df, mut_class==class_2)$A107_score)$p.value
  p_value_A214 <- t.test(filter(df, mut_class==class_1)$A214_score, filter(df, mut_class==class_2)$A214_score)$p.value
  p_value_A218 <- t.test(filter(df, mut_class==class_1)$A218_score, filter(df, mut_class==class_2)$A218_score)$p.value
  print (paste("A107: p-value of diff between", class_1, 'vs', class_2, ':', p_value_A107))
  print (paste("A214: p-value of diff between", class_1, 'vs', class_2, ':', p_value_A214))
  print (paste("A218: p-value of diff between", class_1, 'vs', class_2, ':', p_value_A218))
  }

box_whisk <- function(df, graphname, title){
  df <- df[df$mut_class != 'WT', ]
  textsize <- 7
  p <- ggplot(data=df, aes(x=mut_class, y=param)) +
  geom_boxplot(outlier.shape=NA) +
  geom_sina(alpha=0.6, pch=16, maxwidth=1.2, size=0.2) +
  ggtitle(title) +
  theme_cowplot(12) +
  theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=textsize,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_text(size=textsize,face="bold"),
        axis.title.y=element_text(size=textsize,face="bold")) +
  ylab("Binding score") +
  xlab("") 
  ggsave(graphname, p, width=2, height=2, bg='white', dpi=600)
  }
  
df <- read_tsv('result/S2HR1_scores_common.tsv') %>%
        filter(avg_freq > 0.000015)
print (nrow(df))
t_test(df, 'silent', 'nonsense')
t_test(df, 'silent', 'missense')
t_test(df, 'missense', 'nonsense')
df_plot <- mutate(df, param = A107_score)
box_whisk(df_plot, 'graph/QC_COVA1-07.png','COVA1-07')
df_plot <- mutate(df, param = A214_score)
box_whisk(df_plot, 'graph/QC_COVA2-14.png','COVA2-14')
df_plot <- mutate(df, param = A218_score)
box_whisk(df_plot, 'graph/QC_COVA2-18.png','COVA2-18')
