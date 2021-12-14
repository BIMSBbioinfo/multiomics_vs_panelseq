#!/opt/R/4.0/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("No arguments provided. Aborting", call. = FALSE)
}

# -------------------------------------------------------------------------------------------------
# {SETUP}
# arguments
# 1 - path to working directory
# 2 - intersection toggle
# 3 - filter toggle
# 4 - sum. mutations per gene toggle
# 5 - number of forks
# 6 - Custom gene set to use, if present
# 7 - data specific run code, if present

# paths
p.wd <- args[1]    # /local/abarano/Projects/DrugResponse
p.d <- file.path(p.wd, "data")
if (!dir.exists(file.path(p.wd, "Results"))) {
  dir.create(file.path(p.wd, "Results"))
  dir.create(file.path(p.wd, "Results", "dumpster"))
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
cores <- as.integer(args[5])
# path to gene set provided?
p.gs <- ifelse(args[6] == "NA", NA, args[6])
# run code
run.code <- ifelse(length(args) < 7, paste0(sample(1:1e2, 1), sample(letters, 1)), args[7])



# Report settings
message(paste0("Working directory: ", args[1]))
message(paste0("Intersected features: ", args[2]))
message(paste0("Filtered data: ", args[3]))
message(paste0("Mutatations summarized per gene: ", args[4]))
message(paste0("Number of forked processes: ", args[5]))
message(paste0("Gene set to subset features: ", args[6]))
message(paste0("Run code: ", run.code))


# :°)
set.seed(42)
# load libraries
library(doParallel)
library(foreach)
library(caret)
library(caretEnsemble)
library(tidyverse)
# start the parallelization
Mycluster <- parallel::makeForkCluster(cores)
doParallel::registerDoParallel(Mycluster)


# -------------------------------------------------------------------------------------------------
# {MAIN}
# read the drug response data
drugs <- purrr::map(l.d,
                    ~ data.table::fread(file.path(p.d, "prepared", .x, "drug_response/dr_targeted.tsv")) %>% split(., .$column_name)
)
names(drugs) <- l.d
# -------------------------------------------------------------------------------------------------
# read the data
message("Looking for data...")
if (file.exists(file.path(p.wd, "Results", paste0("PCAembeddings_", run.code, ".RDS")))) {
  message("PCA embddings detected. Reading in...")
  data.pca <- readRDS(file.path(p.wd, "Results", paste0("PCAembeddings_", run.code, ".RDS")))
} else {
  message("PCA embddings not found. Reading in raw data to compute PCA...")
  data <- list()
  data.pca <- list()
  data.loadings <- list()
  for (source in l.d) {
    for (type in d.type) {
      data[[source]][[type]] <- data.table::fread(file.path(p.d, 
                                                            "generated",
                                                            paste0(source, "_", type, 
                                                                   ifelse(toggles$filterData, "_filtered", ""), 
                                                                   ifelse(toggles$intersectFeatures, "_intersected", ""),
                                                                   ifelse(toggles$sumMutPerGene, "_mpg", ""),
                                                                   ".tsv")
                                                            ), 
                                                  data.table = F)
      if (!is.na(p.gs)) {
        if (file.exists(p.gs)) {
          message("Filtering data with the provided gene list...")
          t.gs <- data.table::fread(p.gs, data.table = F)
          v.gs <- t.gs[["gene_id"]]
          v.dnms <- colnames(data[[source]][[type]])
          # filter features based on the provided gene list
          message("Original number of features: ", dim(data[[source]][[type]])[2] - 1)
          data[[source]][[type]] <- data[[source]][[type]][, c(1, which(str_detect(v.dnms, paste0(v.gs, collapse = "|"))))]
          message("Number of features after filtering: ", dim(data[[source]][[type]])[2] - 1)
          
        } else {
          message("No gene list file found. Proceeding with all features...")
        }
      }
      prcomp.obj <- prcomp(data[[source]][[type]][, -1])
      data.pca[[source]][[type]] <- as.data.frame(prcomp.obj$x)
      data.pca[[source]][[type]]$sample_id <- data[[source]][[type]]$sample_id
      data.pca[[source]][[type]] <- data.pca[[source]][[type]][, c("sample_id", paste0("PC", 1:50))]
      # include loadings
      data.loadings[[source]][[type]] <- as.data.frame(prcomp.obj$rotation)
    }
  }
  saveRDS(data.pca, file.path(p.wd, "Results", paste0("PCAembeddings_", run.code, ".RDS")))
  saveRDS(data.loadings, file.path(p.wd, "Results", paste0("PCAloadings_", run.code, ".RDS")))
}


# -------------------------------------------------------------------------------------------------
# prepare training&testing sets
message("Preparing training and testing partitions...")
if (file.exists(file.path(p.wd, "Results", paste0("trainSets_", run.code, ".RDS")))) {
  l.tr <- readRDS(file.path(p.wd, "Results", paste0("trainSets_", run.code, ".RDS")))
  l.tst <- readRDS(file.path(p.wd, "Results", paste0("testSets_", run.code, ".RDS")))
} else {
  l.tr <- list()
  l.tst <- list()
  for(source in l.d) {
    if (source == "PDX_all") {
      drugs.sub <- drugs[[source]]
      drugs.sub.names <- names(drugs.sub)
      for(type in d.type) {
        for(drug in drugs.sub.names) {
          tab.drug <- drugs.sub[[drug]]
          tmp <- inner_join(tab.drug[, c("sample_id", "value")], data.pca[[source]][[type]], by = "sample_id") %>% as.data.frame()
          if(dim(tmp)[1] < 30) {
            next
          }
          l.tr[[source]][[type]][[drug]] <- tmp
        }
      }
    }
    drugs.sub <- drugs[[source]]
    drugs.sub.names <- names(drugs.sub)
    for(type in d.type) {
      for(drug in drugs.sub.names) {
        tab.drug <- drugs.sub[[drug]]
        tmp <- inner_join(tab.drug[, c("sample_id", "value")], data.pca[[source]][[type]], by = "sample_id") %>% as.data.frame()
        if(dim(tmp)[1] < 50) {
          next
        }
        v.tr <- sample(tmp$sample_id, dim(tmp)[1] * 0.8, replace = FALSE)
        v.tst <- tmp$sample_id[!tmp$sample_id %in% v.tr]
        l.tr[[source]][[type]][[drug]] <- tmp[tmp$sample_id %in% v.tr, ]
        l.tst[[source]][[type]][[drug]] <- tmp[tmp$sample_id %in% v.tst, ]
      }
    }
  }
  saveRDS(l.tr, file.path(p.wd, "Results", paste0("trainSets_", run.code, ".RDS")))
  saveRDS(l.tst, file.path(p.wd, "Results", paste0("testSets_", run.code, ".RDS")))
}


# -------------------------------------------------------------------------------------------------
# model specifics
tune_list <- list(
  bridge = caretModelSpec(method = "bridge"),
  blasso = caretModelSpec(method = "blasso",
                          tuneGrid = expand.grid(sparsity = seq(0.1, 0.9, by = 0.15))),
  svmLinear = caretModelSpec(method = "svmLinear2",
                             tuneGrid = expand.grid(cost = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25)),
                             importance = TRUE),
  ranger = caretModelSpec(method = "ranger",
                          tuneGrid = expand.grid(mtry = seq(2, 10, by = 1),
                                                 min.node.size = c(3, 5, 10),
                                                 splitrule = c("variance")),
                          importance = "impurity", 
                          num.threads = 1),
  glmnet = caretModelSpec(method = "glmnet",
                          tuneGrid = expand.grid(alpha = seq(0, 1, length = 30),
                                                 lambda = seq(0.0001, 1, length = 100)))
)
message("Fitting models...")
l.res <- list()
for(source in l.d) {
  for (dataset in d.type) {
    drugs.to.iter <- names(l.tr[[source]][[dataset]])
    if (source == "PDX_all") {
      foreach(drug = drugs.to.iter) %dopar% {
        ## prepare partitions for LOOCV
        subs <- unique(l.tr[[source]][[dataset]][[drug]]$sample_id)
        parts <- vector(mode = "list", length = length(subs))
        for(i in seq_along(subs))
          parts[[i]] <- which(l.tr[[source]][[dataset]][[drug]]$sample_id != subs[i])
        names(parts) <- paste0("LOOCV.", subs)
        ##
        caretEnsemble::caretList(x = l.tr[[source]][[dataset]][[drug]][, -c(1, 2)],
                                 y = l.tr[[source]][[dataset]][[drug]]$value,
                                 trControl = trainControl(
                                   method = "LOOCV",
                                   index = parts,
                                   allowParallel = FALSE,
                                   savePredictions = "final",
                                   returnResamp = "final",
                                   verbose = TRUE
                                 ),
                                 continue_on_fail = FALSE,
                                 tuneList = tune_list
        )
      } -> model.out
      names(model.out) <- drugs.to.iter
    } else {
      foreach(drug = drugs.to.iter) %dopar% {
        caretEnsemble::caretList(x = l.tr[[source]][[dataset]][[drug]][, -c(1, 2)],
                                 y = l.tr[[source]][[dataset]][[drug]]$value,
                                 trControl = trainControl(
                                   method = "cv",
                                   number = 10,
                                   allowParallel = FALSE,
                                   savePredictions = "final",
                                   index = createResample(l.tr[[source]][[dataset]][[drug]]$value, 10),
                                   returnResamp = "final",
                                   verbose = TRUE
                                 ),
                                 continue_on_fail = FALSE,
                                 tuneList = tune_list
        )
      } -> model.out
      names(model.out) <- drugs.to.iter
    }
    #tmp
    l.res[[source]][[dataset]] <- model.out
    saveRDS(l.res, file.path(p.wd, "Results", "dumpster",  paste0("TMP_", source, "_", dataset, "_fitTargeted_", run.code, ".RDS")))
  }
}
saveRDS(l.res, file.path(p.wd, "Results",  paste0("fitTargeted_", run.code, ".RDS")))
parallel::stopCluster(Mycluster)
message("Done!")


