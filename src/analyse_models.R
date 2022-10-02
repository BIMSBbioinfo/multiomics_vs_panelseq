# analyse models built for drug response
## libraries & functions
if (!require("tidyverse")) install.packages("tidyverse") else library(tidyverse)
if (!require("foreach")) install.packages("foreach") else library(foreach)
if (!require("tictoc")) install.packages("tictoc") else library(tictoc)
if (!require("data.table")) install.packages("data.table") else library(data.table)
if (!require("furrr")) install.packages("furrr") else library(furrr)
if (!require("purrr")) install.packages("purrr") else library(purrr)
#utility_functions <- file.path("F_auxiliary.R")
source("src/F_auxiliary.R")

# start the parallelization
message(date()," => assembling model stats")
p.out <- '/data/local/buyar/arcas/panel_sequencing_paper/multiomics_vs_panelseq/test'
outdir <- 'CCLE'

files <- dir(file.path(p.out, outdir), ".caret.RDS$", full.names = T)

cl <- parallel::makeCluster(20)
parallel::clusterExport(cl = cl, varlist = c('files', 'p.out', 'outdir'))
## assemble model stats
dt <- do.call(
  rbind,
  pbapply::pblapply(cl = cl, files,
    function(f) {
      x <- readRDS(f)
      drug <- gsub(".caret.RDS$", "", basename(f))
      dt <- rbind(
        x$panel$panel.stats,
        x$panel$panel.stats.pca,
        x$panel$panel.stats.rf,
        x$panel$panel.stats.rf.pca,
        x$mo$mo.stats,
        x$mo$mo.stats.pca,
        x$mo$mo.stats.rf,
        x$mo$mo.stats.rf.pca
      )
      dt$total_sample_count <- rep(
        c(
          x$mo$total_sample_count,
          x$panel$total_sample_count
        ),
        4
      )
      dt$training_sample_count <- rep(
        c(
          x$mo$training_sample_count,
          x$panel$training_sample_count
        ),
        4
      )
      dt$testing_sample_count <- rep(
        c(
          x$mo$testing_sample_count,
          x$panel$testing_sample_count
        ),
        4
      )
      dt$type <- rep(c("panel", "mo"), each = 4)
      dt$pp <- rep(c("scale+nzv", "scale+nzv+pca"), 4)
      dt$model <- rep(c(rep("glm", 2),rep("rf", 2)),2)
      dt$drug <- rep(drug, 8)
      return(dt)
    }
  )
)
parallel::stopCluster(cl)
saveRDS(dt, file = file.path(p.out, outdir, "caret.stats.RDS"))

## Estimate feature importance glmnet
plan(multisession, workers = 4)
files[1:4] %>%
  future_map(
    .,
    ~ caret::varImp(readRDS(.x)$mo$mo.fit)$importance
  ) -> l.varimp
names(l.varimp) <- gsub(".caret.RDS$", "", basename(files))
                        
purrr::map2(
  l.varimp,
  names(l.varimp),
  ~ .x %>%
    as_tibble(rownames = "Feature") %>%
    mutate(drug = .y)
) %>%
  purrr::reduce(., rbind) -> t.varimp
write_tsv(t.varimp, file.path(p.out, outdir, "varImp.tsv"))


## Estimate feature importance rf
list.files(file.path(p.out, outdir),
           pattern = ".caret.RDS",
           full.names = T
) %>%
  future_map(
    .,
    ~ caret::varImp(readRDS(.x)$mo$mo.fit.rf)$importance
  ) -> l.varimp
names(l.varimp) <- list.files(file.path(p.out, outdir),
                              pattern = ".caret.RDS",
                              full.names = F
) %>%
  basename() %>%
  str_remove(., ".caret.RDS")
purrr::map2(
  l.varimp,
  names(l.varimp),
  ~ .x %>%
    as_tibble(rownames = "Feature") %>%
    mutate(drug = .y)
) %>%
  purrr::reduce(., rbind) -> t.varimp
write_tsv(t.varimp, file.path(p.out, outdir, "varImp_rf.tsv"))