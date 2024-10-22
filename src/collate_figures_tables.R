# Plot manuscript figures from analysis results of CCLE/PDX/beatAML
args <- commandArgs(trailingOnly = T) 

folder <- args[1] #path to folder which contains result folders for CCLE/PDX/beatAML
folder.data <- args[2] #path to a folder that contains all data (including auxiliary, prepared, etc.)

# Check for missing packages, install if needed
options(repos = structure(c(CRAN = "https://ftp.fau.de/cran/")))   # modify according to personal preferences
list_of_pkgs <- c("openxlsx", "ggplot2", "ggpubr", "ggrepel", "ggridges", "patchwork",
                  "data.table", "janitor", "magrittr", "stringr", "dplyr", "tidyr", "grid")
install.packages(list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())])

library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggridges)
library(patchwork)
library(data.table)
library(janitor)
library(magrittr)
library(stringr)
library(dplyr)
library(tidyr)
library(grid)

message(date(), " => collecting results files")

# collate all performance metrics into one table
stats <- do.call(rbind, sapply(simplify = F, c('CCLE', 'PDX', 'beatAML'), function(x) {
  dt <- readRDS(file.path(folder, x, 'caret.stats.RDS'))
  dt$dataset <- x
  return(dt)
}))
stats <- stats %>% dplyr::rename(Rsquared = Rsquare)


# import variable importance metrics 
IMP <- sapply(simplify = F, c('CCLE', 'PDX', 'beatAML'), function(dataset) {
  sapply(simplify = F, c('rf', 'glm', 'svm'), function(a) {
    dt <- data.table::fread(file.path(folder, 
                                      dataset,
                                      paste0('feature_importance.', a, '.tsv')))
  })
})

# import drug classes info
dr.cls <- data.table::fread(file.path(folder.data, "auxiliary", "Repurposing_Hub_export.txt"))
dr.cls$Name <- str_to_upper(dr.cls$Name)
dr.cls <- dr.cls[, .(Name, MOA)]
colnames(dr.cls) <- c("drugName", "MOA")


message(date(), " => making plots")
ggplot2::theme_set(ggpubr::theme_pubclean())
# make a scatter plot of comparing panel vs multiomics, along with a boxplot
# that demonstrates comparison stats
# main figure 1: 
# ccle/beataml/pdx (with RF without pca) (scatter/boxplots)
dt <- stats[model == 'ranger'][ppOpts == 'center+scale+nzv'][dataset %in% c('beatAML', 'CCLE')]
dtc <- dcast.data.table(dt, drugName + dataset ~ type, value.var = 'Rsquared')
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
  labs(color = "Improvement\n(Rsquared)", x = "Panel", y = "Multiomics") +
  facet_wrap(~ dataset, nrow = 1) +
  theme(legend.position = 'right')

# make a figure computing percentage of gex features on top important features for CCLE 
# use "ranger" importance metrics 
# get a summary of top gex 
dt <- stats[model == 'ranger'][ppOpts == 'center+scale+nzv'][dataset == 'CCLE']
dtc <- dcast.data.table(dt, drugName ~ type, value.var = 'Rsquared')
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
    fill = "Rsquared"
  ) +
  theme(legend.position = 'right')


## make PDX plot (boxplot per drug with significance labels)
# plot PDX results for each method 
pdx_plots <- sapply(simplify = F, unique(stats$model), function(m){
  dt <- stats[dataset == 'PDX'][model == m][ppOpts == 'center+scale+nzv']
  dtc <- dcast(dt, drugName + run ~ type, value.var = 'Rsquared')
  stars <- do.call(rbind, lapply(split(dtc, dtc$drugName), function(x) {
    pval <- wilcox.test(x$multiomics, x$panel, alternative = 'greater', 
                        exact = FALSE)[['p.value']]
    data.table('drugName' = x$drugName[1], 
               'pval' = pval,
               'star' = gtools::stars.pval(pval))
  }))
  ggplot(dt, aes( x = drugName, y = Rsquared)) + 
    geom_boxplot(aes(fill = type), outlier.shape = NA) + 
    geom_text(data = stars, aes(x = drugName, y = 0.4, label = star)) +
    scale_fill_brewer(type = 'qual', palette = 6) +
    facet_grid(~ model) +
    labs(x = "Drugs", 
         y = "Rsquared")+
    theme(legend.title = element_blank(), legend.position = 'right') +
    coord_flip()
})

# use "ranger" results as main figure panel for PDX 
p3 <- pdx_plots$ranger

# supp. figure 1: => we say why we chose RF for main figure
# comparison of methods (RF/glmnet/svm) and ppopts (with/without pca)
plots <- lapply(c('center+scale+nzv', 'center+scale+nzv+pca'), function(ds) {
  stats %>% dplyr::filter(dataset %in% c('CCLE', 'beatAML'))%>%
    dplyr::filter(ppOpts == ds)%>%
    ggboxplot(x = 'type', y = 'Rsquared', 
              add = 'jitter', color = 'type') +
    facet_grid(model ~ dataset) + 
    stat_compare_means(paired = T, label.y = 0.4,
                       method.args = list('alternative' = 'greater')) + 
    scale_color_brewer(type = 'qual', palette = 6) +
    theme_bw()+
    theme(legend.position = 'none', axis.title.x = element_blank()) +
    labs(title = paste("Method:",ds))
})

# print 
p <- plots[[1]] + plots[[2]] + pdx_plots$glmnet + pdx_plots$svmRadial + 
  plot_annotation(tag_levels = 'A')
ggsave(filename = file.path(folder, 'figure_S1.pdf'), plot = p, 
       width = 30, height = 15)

# drug classes by improvement 
# assemble a table
dt.dcl <- stats[dataset == "CCLE" & ppOpts == "center+scale+nzv" & model == "ranger"]
dt.dcl <- dt.dcl[, .(Rsquared, drugName, type)]
dt.dcl <- dcast(dt.dcl, drugName ~ type, value.var = "Rsquared")
dt.dcl[, impr := multiomics - panel]
dt.dcl <- dr.cls[dt.dcl, on = "drugName"]
# select 8 drug classes that are the most represented in the current dataset
dt.dcl %>% 
  tabyl(MOA) %>% 
  drop_na() %>%
  arrange(desc(n)) %>% 
  filter(n > 5) %>% 
  pull(MOA) -> v.tmp
# filter the plotting table to include only those
dt.dcl %>% 
  filter(MOA %in% v.tmp) -> t.tmp
# order drug classes based on average improvement in MOs
moa.ord <- t.tmp %>% 
  group_by(MOA) %>% 
  summarise(med_impr = median(impr), 
            avg_impr = mean(impr), .groups = "drop") %>% 
  arrange(avg_impr) %>% 
  pull(MOA)

# assemblt the plot
p4 <- t.tmp %>% 
  mutate(MOA = factor(MOA, levels = c(moa.ord))) %>% 
  ggplot(data = ., aes(x = impr, y = MOA, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 2, 
                               rel_min_height = 0.01, 
                               gradient_lwd = 1,
                               bandwidth = 0.03) + 
  scale_fill_gradient2(low = "black", mid = "gray", high = "red") +
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = 0.053, 
             linetype = "dashed") + 
  scale_x_continuous(expand = c(0.07, -0.01)) +
  labs(x = "Mean Improvement (CCLE,Rsquared)", 
       y = "Drug mechanism of action") +
  guides(fill = guide_colourbar(barwidth = 26, 
                                barheight = 1, 
                                direction = "horizontal", 
                                title = "", label = F, ticks = T)) + 
  theme(legend.position = "bottom", 
        axis.title.x = element_text(vjust = -1.2))

# collate Figure 1 
layout <- "AABB
   AABB 
   CCDD
"

p <- p1 + p2 + p3 + p4 + plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
ggsave(filename = file.path(folder, 'figure_1.pdf'), 
       plot = p, width = 22,
       height = 10)

# supplementary figure 2: 
# drug classes by improvement (glmnet and svm)
lapply(c("glmnet", "svmRadial"), 
       function(model_name) {
         dt.dcl <- stats[dataset == "CCLE" & ppOpts == "center+scale+nzv" & model == model_name]
         dt.dcl <- dt.dcl[, .(Rsquared, drugName, type)]
         dt.dcl <- dcast(dt.dcl, drugName ~ type, value.var = "Rsquared")
         dt.dcl[, impr := multiomics - panel]
         dt.dcl <- dr.cls[dt.dcl, on = "drugName"]
         # select 8 drug classes that are the most represented in the current dataset
         dt.dcl %>% 
           tabyl(MOA) %>% 
           drop_na() %>%
           arrange(desc(n)) %>% 
           filter(n > 5) %>% 
           pull(MOA) -> v.tmp
         # filter the plotting table to include only those
         dt.dcl %>% 
           filter(MOA %in% v.tmp) -> t.tmp
         # order drug classes based on average improvement in MOs
         moa.ord <- t.tmp %>% 
           group_by(MOA) %>% 
           summarise(med_impr = median(impr), 
                     avg_impr = mean(impr), .groups = "drop") %>% 
           arrange(avg_impr) %>% 
           pull(MOA)
         # assemblt the plot
         t.tmp %>% 
           mutate(MOA = factor(MOA, levels = c(moa.ord))) %>% 
           ggplot(data = ., aes(x = impr, y = MOA, fill = stat(x))) + 
           geom_density_ridges_gradient(scale = 2, 
                                        rel_min_height = 0.01, 
                                        gradient_lwd = 1,
                                        bandwidth = 0.03) + 
           scale_fill_gradient2(low = "black", mid = "gray", high = "red") +
           geom_vline(xintercept = 0) + 
           geom_vline(xintercept = 0.053, 
                      linetype = "dashed") + 
           scale_x_continuous(expand = c(0.07, -0.01)) +
           labs(x = paste0("Multiomics improvement (CCLE,", model_name, ",Rsquared)"), 
                y = "Drug mechanism of action") +
           guides(fill = guide_colourbar(barwidth = 26, 
                                         barheight = 1, 
                                         direction = "horizontal", 
                                         title = "", label = F, ticks = T)) + 
           theme(legend.position = "bottom", 
                 axis.title.x = element_text(vjust = -1.2))
       }) -> l.p


# Most important features for venotoclax
# use beatAML data to demonstrate this as it is a AML-drug and shows the top improvement for beatAML but not CCLE.
arrange(IMP$beatAML$rf[IMP$beatAML$rf$Venetoclax > 0, c("feature", "Venetoclax")], desc(Venetoclax)) -> t.vnt.sort
t.vnt.sort <- t.vnt.sort[1:20,]
t.vnt.sort$feature <- str_remove(t.vnt.sort$feature, ".panel")
str_sub(t.vnt.sort$feature, -4, -4) <- "-"
t.vnt.sort$colvar <- str_split_fixed(t.vnt.sort$feature, "-", 2)[, 2]
t.vnt.sort$feature <- str_remove(str_split_fixed(t.vnt.sort$feature, "-", 2)[, 1], "HALLMARK_")
t.vnt.sort$feature <- factor(t.vnt.sort$feature, levels = c(t.vnt.sort$feature))

p.ven <- ggplot(data = t.vnt.sort, aes(x = Venetoclax, y = feature)) + 
  geom_col(aes(fill = colvar), size = 0.5) + 
  scale_fill_brewer(type = 'qual', palette = 6) +
  labs(y = "", x = "Venetoclax feature importance", fill = "Feature type") + 
  theme(legend.position = c(0.75, 0.55))

# collate supplementary figure 2 
layout <- 
  "AAAACC
 BBBBCC"
p <- l.p[[1]] + l.p[[2]] + p.ven + plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
ggsave(filename = file.path(folder, 'figure_S2.pdf'), plot = p, 
       width = 15, height = 8)


message(date(), " => printing supplementary tables to xlsx files")

# supplementary table (export excel)
# 1. performance metrics: ccle/beataml/pdx for rf/glmnet/svm
# 2. variable importance metrics: ccle/beataml/pdx for rf/glmnet/svm

# supp. table 1
# performance metrics
OUT <- createWorkbook()
lapply(unique(stats$dataset), function(Dataset) {
  lapply(unique(stats[dataset == Dataset]$model), function(Model) {
    sname<-paste(Dataset, Model, sep = "_")
    addWorksheet(OUT, sname)
    writeData(OUT, sheet = sname, x = stats[dataset == Dataset][model == Model])
  })
})
saveWorkbook(OUT, file.path(folder, "SupplementaryTable1.xlsx"), overwrite = T)

# supp. table 2
# variable importance metrics
OUT <- createWorkbook( )
lapply(names(IMP), function(dataset) {
  lapply(names(IMP[[dataset]]), function(model) {
    sname<-paste(dataset, "feature_importance", model, sep = "_")
    addWorksheet(OUT, sname)
    writeData(OUT, sheet = sname, x = IMP[[dataset]][[model]])
  })
})
saveWorkbook(OUT, file.path(folder, "SupplementaryTable2.xlsx"), overwrite = T)

# supp. table 3
# improvement across drug classes, per dataset, per model; PDX is omitted because only two drugs were annotated
OUT <- createWorkbook()
lapply(c("CCLE", "beatAML"), function(Dataset) {
  lapply(c("ranger", "glmnet", "svmRadial"), function(Model) {
    dt.dcl <- stats[dataset == Dataset & ppOpts == "center+scale+nzv" & model == Model]
    dt.dcl <- dt.dcl[, .(Rsquared, drugName, type)]
    dt.dcl <- dcast(dt.dcl, drugName ~ type, value.var = "Rsquared")
    dt.dcl[, impr := multiomics - panel]
    dt.dcl[, drugName := str_to_upper(drugName)]
    dt.dcl <- dr.cls[dt.dcl, on = "drugName"]
    dt.dcl <- drop_na(dt.dcl)
    sname <- paste(Dataset, "drug_improv", Model, sep = "_")
    addWorksheet(OUT, sname)
    writeData(OUT, sheet = sname, x = dt.dcl)
  })
})
saveWorkbook(OUT, file.path(folder, "SupplementaryTable3.xlsx"), overwrite = T)


message(date(), " => Finished collating figures and tables")

