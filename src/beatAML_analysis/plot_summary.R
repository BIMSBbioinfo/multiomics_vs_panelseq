# use BeatAML stats to plot summary figures
# make summary figure 

library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)

args <- commandArgs(trailingOnly = T)
stats_file <- args[1]

results <- data.table::fread(stats_file)
results$improvement <- results$Rsquare_multiomics - results$Rsquare_panel

top_drugs <- results[order(improvement, decreasing = T)][1:5]$drug

p1 <- ggplot(results,
             aes(x = Rsquare_panel, y = Rsquare_multiomics)) +
  geom_point(aes(color = improvement), size = 3) + geom_abline(slope = 1) + 
  geom_text_repel(aes(label = ifelse(drug %in% top_drugs, drug, ""))) + 
  coord_fixed() + 
  lims(x = c(0, 0.5), y = c(0, 0.5)) + 
  theme_bw(base_size = 12) +
  scale_color_gradient2(low = 'black', mid = 'gray', high = 'red') +
  labs(color = "Multiomics\nimprovement", x = 'panel', y = 'multiomics') 

dt <- melt.data.table(results[,c('drug', 'Rsquare_multiomics', 'Rsquare_panel')])
dt$variable <- gsub("Rsquare_", "", dt$variable)
dt <- dt[order(variable, decreasing = T)]
p2 <- ggboxplot(dt, x = 'variable', y = 'value', add = 'jitter') + 
  stat_compare_means(paired = T, method.args = list('alternative' = 'greater'),
                     label.x = 2.2, label.y = 0.3)+
  theme_bw(base_size = 12) + 
  theme(axis.title.y = element_blank()) + 
  coord_flip()

p <- cowplot::plot_grid(p1, p2, 
                        ncol = 1, rel_heights = c(3, 1))
ggsave(filename = 'beatAML.plot.pdf', plot = p, width = 5, height = 7) 
