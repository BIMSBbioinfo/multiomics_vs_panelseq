# Plot manuscript figures from analysis results of CCLE/PDX/beatAML
args <- commandArgs(trailingOnly = T) 

folder <- args[1] #path to folder which contains result folders for CCLE/PDX/beatAML

library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(data.table)

message(date(), " => collecting results files")
ggplot2::theme_set(ggpubr::theme_pubclean())

# collate all performance metrics into one table
stats <- do.call(rbind, sapply(simplify = F, c('CCLE', 'PDX', 'beatAML'), function(x) {
  dt <- readRDS(file.path(folder, x, 'caret.stats.RDS'))
  dt$dataset <- x
  return(dt)
}))

# import variable importance metrics 
IMP <- sapply(simplify = F, c('CCLE', 'PDX', 'beatAML'), function(dataset) {
  sapply(simplify = F, c('rf', 'glm', 'svm'), function(a) {
    dt <- data.table::fread(file.path(folder, 
                                      dataset,
                                      paste0('feature_importance.', a, '.tsv')))
  })
})

message(date(), " => making plots")

# make a scatter plot of comparing panel vs multiomics, along with a boxplot
# that demonstrates comparison stats
# main figure 1: 
# ccle/beataml/pdx (with RF without pca) (scatter/boxplots)
dt <- stats[model == 'ranger'][ppOpts == 'center+scale+nzv'][dataset %in% c('beatAML', 'CCLE')]
dtc <- dcast.data.table(dt, drugName + dataset ~ type, value.var = 'Rsquare')
# get top 5 drugs per dataset
top_drugs <- do.call(rbind, lapply(split(dtc, dtc$dataset), function(x) {
  x[order(multiomics - panel, decreasing = T),][1:5]
}))


p1 <- ggscatter(dtc, x = 'panel', y = 'multiomics') +
  geom_point(aes(color = multiomics - panel), size = 2) +
  geom_abline(slope = 1) +
  lims(x = c(0, 0.6), y = c(0, 0.6)) +
  coord_fixed() + 
  geom_text_repel(data = top_drugs, aes(label = drugName), size = 3) + 
  scale_color_gradient2(low = "black", mid = "gray", high = "red") +
  labs(color = "Improvement", x = "panel", y = "multiomics") +
  facet_wrap(~ dataset, nrow = 1) +
  theme(legend.position = 'right')

# make a figure computing percentage of gex features on top important features for CCLE 
# use "ranger" importance metrics 
# get a summary of top gex 
dt <- stats[model == 'ranger'][ppOpts == 'center+scale+nzv'][dataset == 'CCLE']
dtc <- dcast.data.table(dt, drugName ~ type, value.var = 'Rsquare')
dtc$improvement <- dtc$multiomics - dtc$panel
# compute percentage of "gex" features in top informative features for CCLE
# use ranger importance metrics 
I <- IMP$CCLE$rf
features <- do.call(rbind, lapply(colnames(I[,-1]), function(drug) {
  # sort features by importance, compute percentage of gex in top 100
  top <- I[,c('feature', drug),with = F][order(I[[drug]], decreasing = T)][1:100]
  data.table('drugName' = drug, 'count_gex' = length(grep('gex', top$feature)), 
             'count_panel' = length(grep('panel', top$feature)))
}))
dtc <- merge.data.table(dtc, features, by = 'drugName')
p2 <- ggscatter(dtc, x = 'improvement', y = 'count_gex', fill = 'multiomics', 
                shape = 21, size = 3, cor.coef = T) + 
  #  geom_point(aes(fill = multiomics), shape = 21, size = 3) +
  geom_vline(
    xintercept = c(0, mean(dtc$improvement, na.rm = T)),
    linetype = c("solid", "dashed"), color = c("black", "black"), size = c(0.5, 1)
  ) +
  geom_text(
    data = data.frame(
      x = mean(dtc$improvement, na.rm = T) + 0.01, y = 1.5,
      text = paste0(
        "Mean improvement = ",
        round(mean(dtc$improvement, na.rm = T), 3)
      )
    ),
    aes(x, y, label = text), hjust = 'left'
  ) +
  scale_fill_gradient2(mid = "gray", high = "red") +
  labs(
    x = "Multiomics improvement (CCLE)",
    y = "% of gene expression features\namong top 100 Features",
    fill = "Rsquare"
  ) +
  theme(legend.position = 'right')



## make PDX plot (boxplot per drug with significance labels)
dt <- stats[dataset == 'PDX'][model == 'ranger'][ppOpts == 'center+scale+nzv']
dtc <- dcast(dt, drugName + run ~ type, value.var = 'Rsquare')
stars <- do.call(rbind, lapply(split(dtc, dtc$drugName), function(x) {
  pval <- wilcox.test(x$multiomics, x$panel, alternative = 'greater')[['p.value']]
  data.table('drugName' = x$drugName[1], 
             'pval' = pval,
             'star' = gtools::stars.pval(pval))
}))
p3 <- ggplot(dt, aes( x = drugName, y = Rsquare)) + 
  geom_boxplot(aes(fill = type), outlier.shape = NA) + 
  geom_text(data = stars, aes(x = drugName, y = 0.5, label = star)) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  facet_grid(~ dataset) +
  theme(legend.title = element_blank(), legend.position = 'top') +
  coord_flip()


# collate Figure 1 
layout <- 
"AAACC
AAACC
BBBCC
"
p <- p1 + p2 + p3 + plot_layout(design = layout) + plot_annotation(tag_levels = 'A')

ggsave(filename = 'figure1.pdf', plot = p, width = 13.76, height = 7.44)

# supp. figure 1: => we say why we chose RF for main figure
# comparison of methods (RF/glmnet/svm) and ppopts (with/without pca)
plots <- lapply(c('CCLE', 'beatAML'), function(ds) {
  ggboxplot(stats[dataset == ds], x = 'type', y = 'Rsquare', 
            add = 'jitter', color = 'type') +
    facet_grid(ppOpts ~ model) + stat_compare_means() + 
    scale_color_brewer(type = 'qual', palette = 6) +
    theme(legend.position = 'none', axis.title.x = element_blank()) + 
    labs(title = ds)
})

p <- plots[[1]] + plots[[2]]
ggsave(filename = 'figure_S1.pdf', plot = p, 
       width = 14.7, height = 7.44)

# main figure 2: 
# drug classes by improvement 



message(date(), " => printing supplementary tables to xlsx files")

# supplementary table (export excel)
# 1. performance metrics: ccle/beataml/pdx for rf/glmnet/svm
# 2. variable importance metrics: ccle/beataml/pdx for rf/glmnet/svm

# supp. table 1
# performance metrics
OUT <- createWorkbook()
lapply(unique(stats$dataset), function(dataset) {
  lapply(unique(stats[dataset == dataset]$model), function(model) {
    sname<-paste(dataset, model, sep = "_")
    addWorksheet(OUT, sname)
    writeData(OUT, sheet = sname, x = stats[dataset == dataset][model == model])
  })
})
saveWorkbook(OUT, "SupplementaryTable1.xlsx")

# supp. table 2
# variable importance metrics
OUT <- createWorkbook()
lapply(names(IMP), function(dataset) {
  lapply(names(IMP[[dataset]]), function(model) {
    sname<-paste(dataset, "feature_importance", model, sep = "_")
    addWorksheet(OUT, sname)
    writeData(OUT, sheet = sname, x = IMP[[dataset]][[model]])
  })
})
saveWorkbook(OUT, "SupplementaryTable2.xlsx")

message(date(), " => Finished collating figures and tables")






