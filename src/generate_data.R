#!/opt/R/4.0/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No arguments provided. Aborting", call. = FALSE)
}

# -------------------------------------------------------------------------------------------------
# {SETUP}
# arguments
# 1 - path to working directory
# 2 - intersection toggle

# paths
p.wd <- args[1]    # /local/abarano/Projects/DrugResponse
p.d <- file.path(p.wd, "data")
if (!dir.exists(file.path(p.d, "generated"))) {
  dir.create(file.path(p.d, "generated"))
}


# data specifics
l.d <- c("CCLE", "PDX_all")
d.type <- c("multiomics", "panelseq", "exomeseq")


# settings
toggles <- list()
toggles$intersectFeatures <- as.logical(args[2])
toggles$filterData <- as.logical(args[3])
toggles$sumMutPerGene <- as.logical(args[4])

# number of the cores to be used
cores <- 6

# :°)
set.seed(42)

# load libraries
library(qsmooth)
library(tidyverse)
library(furrr)
library(FeatureSelection)
library(foreach)
library(doParallel)

source(file.path(p.wd, "src/functions.R"))


# start the parallelization
Mycluster <- parallel::makeForkCluster(cores)
doParallel::registerDoParallel(Mycluster)

# -------------------------------------------------------------------------------------------------
# {MAIN}
  for (data in d.type) {
    results <- list()
    if (data %in% c("multiomics", "exomeseq")) {
      k <- "all_intersect_genes"
    } else {
      k <- "oncoKB_intersect_genes"
    }
    cnv <- list()
    mut <- list()
    gex <- list()
    for (d in l.d) {
      if (!toggles$filterData) {
        cnv[[d]] <- loadCNV(file.path(p.d, "prepared"), d, k)
        mut[[d]] <- loadMut(file.path(p.d, "prepared"), d, k)
      } else {
        cnv[[d]] <- filtCNV(cnv[[d]])
        mut[[d]] <- filtMut(mut[[d]], var_coef = -2, var_intercept = 0.01)
      }
      if (toggles$sumMutPerGene) {
        mut[[d]] <- sumMutPerGene(mut[[d]])
      }
    }
    if (data == "multiomics") {
      for (d in l.d) {
        if (!toggles$filterData) {
          gex[[d]] <- loadGEX(file.path(p.d, "prepared"), d, k)
        } else {
          gex[[d]] <- filtGEX(gex[[d]], var_coef = -2 / 3, var_intercept = 1.5)
        }
      }
      for (d in l.d) {
        if (toggles$intersectFeatures) {
          if (!toggles$filterData) {
            # tmp.resList <- purrr::map(c("cnv", "mut", "gex"), 
            #                           ~ eval(parse(text = .x))$data_raw[])
            tmp.resList <- list(
              cnv[[d]] <- cnv[[d]]$data_raw[, purrr::map(cnv, ~ colnames(.x$data_raw)) %>% purrr::reduce(., intersect)],
              mut[[d]] <- mut[[d]]$data_raw[, purrr::map(cnv, ~ colnames(.x$data_raw)) %>% purrr::reduce(., intersect)],
              gex[[d]] <- gex[[d]]$data_raw[, purrr::map(cnv, ~ colnames(.x$data_raw)) %>% purrr::reduce(., intersect)]
            )           
          } else {
            tmp.resList <- list(
              cnv[[d]] <- cnv[[d]]$data_varFiltered[, purrr::map(cnv, ~ colnames(.x$data_varFiltered)) %>% purrr::reduce(., intersect)],
              mut[[d]] <- mut[[d]]$data_varFiltered[, purrr::map(cnv, ~ colnames(.x$data_varFiltered)) %>% purrr::reduce(., intersect)],
              gex[[d]] <- gex[[d]]$data_varFiltered[, purrr::map(cnv, ~ colnames(.x$data_varFiltered)) %>% purrr::reduce(., intersect)]
            )            
          }

        } else {
          if (!toggles$filterData) {
            tmp.resList <- list(
              cnv[[d]]$data_raw,
              mut[[d]]$data_raw,
              gex[[d]]$data_raw
            )
          } else {
            tmp.resList <- list(
              cnv[[d]]$data_varFiltered,
              mut[[d]]$data_varFiltered,
              gex[[d]]$data_varFiltered
            )
          }
        }
        results[[d]] <- tmp.resList %>%
          future_map2(
            .,
            c("cnv", "mut", "gex"),
            ~ .x %>%
              `colnames<-`(paste0(colnames(.x), "_", .y)) %>%
              as_tibble(., rownames = "sample_id")
          ) %>%
          purrr::reduce(., inner_join, by = "sample_id")
        # save the results
        write_tsv(results[[d]], 
                  file.path(p.d, 
                            "generated",
                            paste0(d, "_", data, 
                                   ifelse(toggles$filterData, "_filtered", ""), 
                                   ifelse(toggles$intersectFeatures, "_intersected", ""),
                                   ifelse(toggles$sumMutPerGene, "_mpg", ""),
                                   ".tsv")
                            )
                  )
      }
      gc()
    } else {
      for (d in l.d) {
        if (toggles$intersectFeatures) {
          if (!toggles$filterData) {
            tmp.resList <- list(
              cnv[[d]] <- cnv[[d]]$data_raw[, purrr::map(cnv, ~ colnames(.x$data_raw)) %>% purrr::reduce(., intersect)],
              mut[[d]] <- mut[[d]]$data_raw[, purrr::map(cnv, ~ colnames(.x$data_raw)) %>% purrr::reduce(., intersect)]
            )  
          } else {
            tmp.resList <- list(
              cnv[[d]] <- cnv[[d]]$data_varFiltered[, purrr::map(cnv, ~ colnames(.x$data_varFiltered)) %>% purrr::reduce(., intersect)],
              mut[[d]] <- mut[[d]]$data_varFiltered[, purrr::map(cnv, ~ colnames(.x$data_varFiltered)) %>% purrr::reduce(., intersect)]
            )  
          }
        } else {
          if (!toggles$filterData) {
            tmp.resList <- list(
              cnv[[d]]$data_raw,
              mut[[d]]$data_raw
            )
          } else {
            tmp.resList <- list(
              cnv[[d]]$data_varFiltered,
              mut[[d]]$data_varFiltered
            )
          }
        }
        results[[d]] <- tmp.resList %>%
          future_map2(
            .,
            c("cnv", "mut"),
            ~ .x %>%
              `colnames<-`(paste0(colnames(.x), "_", .y)) %>%
              as_tibble(., rownames = "sample_id")
          ) %>%
          purrr::reduce(., inner_join, by = "sample_id")
        # save the results
        write_tsv(results[[d]], 
                  file.path(p.d, 
                            "generated",
                            paste0(d, "_", data, 
                                   ifelse(toggles$filterData, "_filtered", ""), 
                                   ifelse(toggles$intersectFeatures, "_intersected", ""),
                                   ifelse(toggles$sumMutPerGene, "_mpg", ""),
                                   ".tsv")
                            )
                  )
      }
      gc()
    }
  }
parallel::stopCluster(Mycluster)