args = commandArgs(trailingOnly = TRUE)
require(tidyverse)
require(janitor)
require(patchwork)


## Setup
p.script.dir <- getwd()
p.parent.dir <- dirname(p.script.dir)
path.in <- as.character(args[1])
path.out <- as.character(ifelse(length(args) < 2, p.parent.dir, args[2]))
## Data - to recreate figures from the results
res <- list()
res[["CCLE"]] <- readRDS(file.path(path.in, "88x_caretRes/caret.stats.RDS"))
res[["CCLE_varImp"]] <- read_tsv(file.path(path.in, "88x_caretRes/varImp.tsv")
res[["PDX.Rep"]] <- readRDS(file.path(path.in, "PDX.Rep_caret.stats.RDS")

## Plotting
### 1
res$CCLE %>% 
  dplyr::select(-RMSE, -COR) %>% 
  dplyr::filter(type == "mo") %>% 
  pivot_wider(, names_from = pp, values_from = Rsquare) -> p.mo
res$CCLE %>% 
  dplyr::select(-RMSE, -COR) %>% 
  dplyr::filter(type == "panel") %>% 
  pivot_wider(, names_from = pp, values_from = Rsquare) -> p.panel
(plot_pcaImp(p.mo) | plot_pcaImp(p.panel)) -> p.1
ggsave(plot = p.1,
       filename = file.path(path.out, "CCLE_Improv.PCA.pdf"), 
       device = "pdf", width = 10, height = 4)
### 2
res$CCLE %>% 
  dplyr::select(-RMSE, -COR) %>% 
  dplyr::filter(pp == "scale+nzv") %>% 
  pivot_wider(, names_from = type, values_from = Rsquare) -> p.npca
res$CCLE %>% 
  dplyr::select(-RMSE, -COR) %>% 
  dplyr::filter(pp == "scale+nzv+pca") %>% 
  pivot_wider(, names_from = type, values_from = Rsquare) -> p.pca
(plot_moImp(p.npca) | plot_moImp(p.pca)) -> p.2
ggsave(plot = p.2,
       filename = file.path(path.out, "CCLE_Improv.MO.pdf"), 
       device = "pdf", width = 10, height = 4)

### 3 Process PDX repeated
res[["PDX.Rep"]] %>% 
  filter(type == "mo", pp == "scale+nzv", !is.na(Rsquare)) %>% 
  group_by(drug) %>% 
  summarize(rsq_med = median(Rsquare), .groups = "drop") %>% 
  arrange(rsq_med) %>% 
  pull(drug) -> dr.ord.pdx
res[["PDX.Rep"]] %>% 
  mutate(drug = factor(drug, levels = dr.ord.pdx)) %>% 
  ggplot(data = ., aes(x = Rsquare, y = drug)) + 
  geom_boxplot(aes(fill = type)) +
  scale_fill_manual(values = alpha(c("orangered", "steelblue"), 0.67)) + 
  lims(x = c(0, 0.6)) +
  theme_bw(base_size = 12) + 
  labs(x = "R-squared", y = "Drug name") + 
  facet_grid(. ~ pp) -> p.3
ggsave(plot = p.3,
       filename = file.path(path.out, "PDX_Improv.MO.pdf"), 
       device = "pdf", width = 6.5, height = 4)

### 4
res[["CCLE_varImp"]] %>% 
  group_by(drug) %>% 
  arrange(desc(Overall)) %>% 
  dplyr::slice(1:100) %>% 
  summarise(cts.gex = sum(str_count(Feature, "gex")),
            cts.mut = sum(str_count(Feature, "mut")),
            cts.cnv = sum(str_count(Feature, "cnv")),
            .groups = "drop") %>% 
  inner_join(p.npca, ., by = "drug") -> t.plot.varImp
ggplot(data = t.plot.varImp,
       aes(x = mo - panel, 
           y = 100 * cts.gex / (cts.gex + cts.mut + cts.cnv))
       ) + 
  geom_point(aes(fill = mo), shape = 21, size = 3) +
  geom_vline(xintercept = c(0, mean(t.plot.varImp$mo - t.plot.varImp$panel, na.rm = T)), 
             linetype = c("solid", "dashed"), color = c("black", "black"), size = c(0.5, 1)) + 
  geom_text(data = tibble(x = mean(t.plot.varImp$mo - t.plot.varImp$panel, na.rm = T) + 0.08, y = 1.5, 
                          text = paste0("Mean improvement = ",
                                        round(mean(t.plot.varImp$mo - t.plot.varImp$panel, na.rm = T), 3))), 
            aes(x, y, label = text)) +
  scale_fill_gradient2(mid = "white", high = "orangered") + 
  #scale_color_gradientn(colours = terrain.colors(100)) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = 0.6, label.y.npc = 0.9) +
  labs(x = "Multi-omics improvement\nR^2 (multi-omics) - R^2 (panel-seq)", 
       y = "% of gene expression features\namong top 100 PCs",
       fill = "R^2 (multi-omics)") +
  theme_bw(base_size = 12) -> p.4
ggsave(plot = p.4,
       filename = file.path(path.out, "CCLE_Improv.Varimp.pdf"), 
       device = "pdf", width = 5.8, height = 4)









