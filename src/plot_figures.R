#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

require(tidyverse)
require(janitor)
require(patchwork)

## auxiiary
plot_pcaImp <- function(data_melted, base.size = 12) {
  (
    (ggplot(data = data_melted, aes(y = `scale+nzv+pca`, x = "")) + 
       geom_boxplot() + 
       lims(y = c(0, 0.5)) +
       theme_bw(base_size = base.size) + labs(x = "")
    ) +
      (ggplot(data = data_melted, aes(x = `scale+nzv`, y = `scale+nzv+pca`)) + 
         geom_point(aes(color = `scale+nzv+pca` - `scale+nzv`)) + 
         geom_abline(intercept = 0, slope = 1) + 
         lims(x = c(0, 0.5), y = c(0, 0.5)) +
         scale_color_gradientn(colours = terrain.colors(100)) +
         theme_bw(base_size = base.size) + labs(x = "", y = "", color = "PCA\nimprovement")
      ) +
      (plot_spacer()
      ) + 
      (ggplot(data = data_melted, aes(x = `scale+nzv`, y = "")) + 
         geom_boxplot() + 
         lims(x = c(0, 0.5)) +
         theme_bw(base_size = base.size) + labs(y = "")
      )
  ) + plot_layout(widths = c(0.2, 1),
                  heights = c(1, 0.2))
}
plot_moImp <- function(data_melted, base.size = 12) {
  (
    (ggplot(data = data_melted, aes(y = mo, x = "")) + 
       geom_boxplot() + 
       lims(y = c(0, 0.5)) +
       theme_bw(base_size = base.size) + labs(x = "")
    ) +
      (ggplot(data = data_melted, aes(x = panel, y = mo)) + 
         geom_point(aes(color = mo - panel)) + 
         geom_abline(intercept = 0, slope = 1) + 
         lims(x = c(0, 0.5), y = c(0, 0.5)) +
         scale_color_gradientn(colours = terrain.colors(100)) +
         theme_bw(base_size = base.size) + labs(x = "", y = "", color = "Multiomics\nimprovement")
      ) +
      (plot_spacer()
      ) + 
      (ggplot(data = data_melted, aes(x = panel, y = "")) + 
         geom_boxplot() + 
         lims(x = c(0, 0.5)) +
         theme_bw(base_size = base.size) + labs(y = "")
      )
  ) + plot_layout(widths = c(0.2, 1),
                  heights = c(1, 0.2))
}


## Setup
path.out <- args[1]


## Data
res <- list()
res[["CCLE"]] <- readRDS("/local/abarano/Projects/DrugResponse/Results/88x_caretRes/caret.stats.RDS")
res[["PDX"]] <- readRDS("/local/abarano/Projects/DrugResponse/Results/81k_caretRes/caret.stats.RDS")


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
ggsave(plot = p.1,
       filename = file.path(path.out, "CCLE_Improv.MO.pdf"), 
       device = "pdf", width = 10, height = 4)


### 3
res$PDX %>% 
  dplyr::select(-RMSE, -COR) %>% 
  dplyr::filter(pp == "scale+nzv") %>% 
  pivot_wider(, names_from = type, values_from = Rsquare) %>% 
  mutate(panel = ifelse(is.na(panel), 0, panel),
         mo = ifelse(is.na(mo), 0, mo)) -> p.pdx
p.pdx.ordr <- p.pdx %>% arrange(mo) %>% pull(drug)
res$PDX %>% 
  dplyr::select(-RMSE, -COR) %>% 
  dplyr::filter(pp == "scale+nzv+pca") %>% 
  pivot_wider(, names_from = type, values_from = Rsquare) %>% 
  mutate(panel = ifelse(is.na(panel), 0, panel),
         mo = ifelse(is.na(mo), 0, mo)) -> p.pdx.pca

(ggplot(data = p.pdx %>% mutate(drug = factor(drug, levels = p.pdx.ordr))) + 
    geom_segment(aes(x = panel, xend = mo, y = drug, yend = drug), color = "grey50") +
    geom_point(aes(x = panel, y = drug), color = "steelblue") + 
    geom_point(aes(x = mo, y = drug,), color = "orangered") + 
    lims(x = c(0, 0.25)) +
    theme_bw() + labs(x = "R-squared", y = "Drug name")) | 
  (ggplot(data = p.pdx.pca %>% mutate(drug = factor(drug, levels = p.pdx.ordr))) + 
     geom_segment(aes(x = panel, xend = mo, y = drug, yend = drug), color = "grey50") +
     geom_point(aes(x = panel, y = drug), color = "steelblue") + 
     geom_point(aes(x = mo, y = drug,), color = "orangered") + 
     lims(x = c(0, 0.25)) +
     theme_bw() + labs(x = "R-squared", y = "Drug name")) -> p.3
ggsave(plot = p.3,
       filename = file.path(path.out, "PDX_Improv.MO.pdf"), 
       device = "pdf", width = 6, height = 4)


### 3
t.varimp <- read_tsv("/local/abarano/Projects/DrugResponse/Results/88x_caretRes/varImp.tsv")
t.varimp %>% 
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









