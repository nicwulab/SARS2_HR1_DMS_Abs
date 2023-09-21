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

plot_score_vs_pos <- function(df, r_table, graphname, xlab, ylab){
   textsize <- 7
   ymax <- 1.2
   colorscale  <- c(brewer.pal(9,"Set1"))[4:5]
   p <-  ggplot() +
     #geom_point(data=df, aes(x=pos, y=mean_adjust_bind), color='gray', alpha=1, pch=16, size=0.3) +
     geom_area(data=df, aes(x=pos, y=mean_adjust_bind), fill='black', alpha=1, color='black') +
     geom_line(data=df, aes(x=pos, y=mean_adjust_bind), color='black', linewidth=0.2) +
     geom_abline(intercept = 0, slope = 0, col='black', linewidth=0.5) +
     geom_rect(data=r_table, mapping=aes(xmin=start, xmax=end, ymin=ymax*0.9, ymax=ymax, fill=region)) + 
     scale_fill_manual(values=colorscale,drop=FALSE) +
     theme_cowplot(12) +
     theme(plot.title=element_text(size=textsize+1,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
           axis.text=element_text(size=textsize,face="bold",colour = 'black'),
           axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
           axis.text.y=element_text(hjust=1,vjust=0.5,colour = 'black'),
           axis.title.x=element_text(size=textsize,face="bold"),
           axis.title.y=element_text(size=textsize,face="bold"),
           legend.position = "none",
           legend.title    = element_text(size=textsize,face="bold"),
           legend.text=element_text(size=textsize,face="bold"),
           legend.justification='center',
           legend.key.size = unit(0.3,"line")) +
     ylab(ylab) +
     xlim(883,1034) +
     xlab(xlab)
   ggsave(graphname,p,width=2.5, height=2, bg='white', dpi=600)
 }

df <- read_tsv('result/S2HR1_adj_score_by_resi.tsv')
r_table <- read_tsv('data/regions.tsv')
plot_score_vs_pos(df, r_table, 'graph/adj_score_vs_pos.png', 'residue', 'mean(adj score)')
