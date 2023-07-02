library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggbreak)
library(patchwork)
require(cowplot)

plot_MFI_change <- function(stats, df, graphname, title){
  palette <- c('grey',qualpal(n = 8, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  palette <- c(palette[8], palette[3], palette[4], palette[6])
  textsize <- 7
  p <- ggplot() +
    geom_bar(stats, mapping = aes(x = x_name, y = mean, colour = color, fill = color),
	     stat = 'identity', width = 0.8, linewidth=0.5, alpha = 0.7) +
    geom_point(df, mapping = aes(x = name, y = value, colour = color, fill = color), size = 0.4, shape = 1, stroke = 0.8) +
    scale_color_manual('', values=palette,drop=FALSE) +
    scale_fill_manual('', values=palette,drop=FALSE) +
    ylim(-0.1, 7) +
    xlab("") +
    ylab("MFI fold change") +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
	  axis.text=element_text(size=textsize,face="bold",colour = 'black'),
	  axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
	  axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
	  axis.title.x=element_blank(),
	  axis.title.y=element_text(size=textsize,face="bold"),
          axis.line.y.right = element_blank(),
	  axis.text.y.right = element_blank(),
	  axis.ticks.y.right = element_blank(),
	  legend.position = "none") +
    ggtitle(title)
  if (graphname!='graph/expr_esc.png'){
    p <- p+scale_y_break(c(1.3, 1.79), scales=0.5, ticklabels=c(2, 4, 6)) 
    }
  if (graphname=='graph/expr_esc.png'){
    ggsave(graphname, width = 1.4, height = 1.82, plot = p, dpi=600, bg='white')
    }
  else{
    ggsave(graphname, width = 1.7, height = 2, plot = p, dpi=600, bg='white')
    }
  }

coloring <- function(mut){
  if (mut %in% c("T961F", "V987C", "Q1010W", "Q1005R")){return ("Low binding")}
  else if (mut %in% c("V915H")){return ("High binding")}
  else if (mut %in% c("WT")){return ("WT")}
  else if (mut %in% c("D950N", "Q954H", "N969K", "L981F", "S982A", "T1027I")){return ("Natural variants")}
  else {return ('Others')}
  }

#print (paste('p-value between WT and Q1010W: ', t.test(df$WT, df$Q1010W, alternative = c('two.sided'))$p.value))
wrapper <- function(infile, graphname, title){
  color_levels <- c('High binding', 'Low binding', 'Natural variants','WT')
  df <- read.csv(file = infile, header = TRUE) 
  avg <- colMeans(df)
  sd <- sapply(df, sd)
  stats <- data.frame(x_name = colnames(df), mean=avg, sd = sd) %>%
             mutate(x_name=factor(x_name, levels=x_name)) %>%
             mutate(color=mapply(coloring, x_name)) %>%
             mutate(color=factor(color, levels=color_levels))
  df <- pivot_longer(df, cols = everything()) %>%
          mutate(color=mapply(coloring, name)) %>%
          mutate(color=factor(color, levels=color_levels))
  plot_MFI_change(stats, df, graphname, title)
  }

wrapper('data/A107_esc.csv','graph/A107_esc.png', 'COVA1-07')
wrapper('data/A214_esc.csv','graph/A214_esc.png', 'COVA2-14')
wrapper('data/A218_esc.csv','graph/A218_esc.png', 'COVA2-18')
wrapper('data/expr_esc.csv','graph/expr_esc.png', 'CC12.3')

