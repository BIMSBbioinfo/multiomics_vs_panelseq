args <- commandArgs(trailingOnly = TRUE)

# {SETUP}
## Paths
p.parent.dir <- getwd()
dset <- as.character(args[1]) ## Either CCLE or PDX
p.data <- ifelse(length(args) >= 2, as.character(args[2]), file.path(p.parent.dir, "data")) ## Path to a folder that contains "prepared" data
p.out <- ifelse(length(args) >= 3, as.character(args[3]), p.parent.dir) ## Path to a folder where the models and their summary statistics will be written
## libraries & functions
if (!require("tidyverse")) install.packages("tidyverse") else library(tidyverse)
if (!require("foreach")) install.packages("foreach") else library(foreach)
if (!require("tictoc")) install.packages("tictoc") else library(tictoc)
if (!require("data.table")) install.packages("data.table") else library(data.table)
if (!require("furrr")) install.packages("furrr") else library(furrr)
if (!require("purrr")) install.packages("purrr") else library(purrr)
#utility_functions <- file.path("F_auxiliary.R")
source("src/F_auxiliary.R")
## read in data
dat <- readRDS(file.path(p.data, "prepared", paste0("data_", dset, ".RDS")))
#dat <- readRDS(file.path(p.data, paste0("data_", dset, ".RDS")))
dr <- data.table::fread(file.path(p.data, "Raw", dset, "drug_response.tsv")) # get drug response data

# {MAIN}
# df: samples on rows, features on columns, includes outcome variable 'y'
train_caret.rf <- function(df, ppOpts = c("center", "scale")) {
  require(caret)
  tgrid <- expand.grid(
    .mtry = seq(from = 10, to = round(sqrt(ncol(df))), 10),
    .splitrule = "variance",
    .min.node.size = c(5, 10, 15)
  )
  model_caret <- train(y ~ .,
    data = df,
    method = "ranger",
    trControl = trainControl(
      method = "repeatedcv", number = 5,
      verboseIter = T, repeats = 3
    ),
    tuneGrid = tgrid,
    preProcess = ppOpts,
    num.trees = 500,
    importance = "permutation"
    #num.threads = 2
  )
  return(model_caret)
}
train_caret.glm <- function(df, ppOpts = c("center", "scale")) {
  require(caret)
  tgrid <- expand.grid(
    .alpha = seq(0, 1, length = 10),
    .lambda = seq(0.0001, 1, length = 20)
  )
  model_caret <- train(y ~ .,
    data = df,
    method = "glmnet",
    trControl = trainControl(
      method = "repeatedcv", number = 5,
      verboseIter = T, repeats = 3
    ),
    tuneGrid = tgrid,
    preProcess = ppOpts
  )
  return(model_caret)
}
run_caret <- function(dat, dr, drugName) {
  selected <- intersect(dr[column_name == drugName]$sample_id, colnames(dat$mut))
  colData <- dr[column_name == drugName][match(selected, sample_id)]

  # change feature names to avoid overlaps from different layers
  dat.raw <- sapply(simplify = F, names(dat), function(x) {
    m <- dat[[x]][, selected]
    rownames(m) <- paste0(rownames(m), "_", x)
    return(m)
  })
  train_samples <- sample(selected, round(length(selected) * 0.7))
  test_samples <- setdiff(selected, train_samples)

  y.train <- colData[match(train_samples, sample_id)]$value
  y.test <- colData[match(test_samples, sample_id)]$value

  # compute results for mut+cnv features (exomeseq) -------------------
  # message(date(), " => processing exomseq")
  # ex.train <- data.frame(do.call(cbind, lapply(dat.raw[c('mut', 'cnv')],
  #                                              function(x) t(x[,train_samples]))),
  #                        check.names = F)
  # ex.train$y <- y.train
  # ex.test <- data.frame(do.call(cbind, lapply(dat.raw[c('mut', 'cnv')],
  #                                             function(x) t(x[,test_samples]))),
  #                       check.names = F)
  # ex.test$y <- y.test
  # ex.fit.rf.pca <- train_caret.rf(ex.train, c("center", "scale", "nzv", "pca"))

  # compute results for panel --------------------------------------------------
  message(date(), " => processing panel")
  panel.train <- data.frame(do.call(
    cbind,
    lapply(
      dat.raw[c("mut.panel", "cnv.panel")],
      function(x) t(x[, train_samples])
    )
  ),
  check.names = F
  )
  panel.train$y <- y.train
  panel.test <- data.frame(do.call(cbind, lapply(
    dat.raw[c("mut.panel", "cnv.panel")],
    function(x) t(x[, test_samples])
  )),
  check.names = F
  )
  panel.test$y <- y.test
  panel.fit.rf <- train_caret.rf(panel.train, c("center", "scale", "nzv"))
  panel.fit.rf.pca <- train_caret.rf(panel.train, c("center", "scale", "nzv", "pca"))
  panel.fit.glm <- train_caret.glm(panel.train, c("center", "scale", "nzv"))
  panel.fit.glm.pca <- train_caret.glm(panel.train, c("center", "scale", "nzv", "pca"))



  # compute results for mut+cnv+gex features (mo=multiomics) -------------------
  message(date(), " => processing multiomics")
  mo.train <- data.frame(do.call(cbind, lapply(
    dat.raw[c("mut.panel", "cnv.panel", "gex")],
    function(x) t(x[, train_samples])
  )),
  check.names = F
  )
  mo.train$y <- y.train
  mo.test <- data.frame(do.call(cbind, lapply(
    dat.raw[c("mut.panel", "cnv.panel", "gex")],
    function(x) t(x[, test_samples])
  )),
  check.names = F
  )
  mo.test$y <- y.test
  mo.fit.rf <- train_caret.rf(mo.train, c("center", "scale", "nzv"))
  mo.fit.rf.pca <- train_caret.rf(mo.train, c("center", "scale", "nzv", "pca"))
  mo.fit.glm <- train_caret.glm(mo.train, c("center", "scale", "nzv"))
  mo.fit.glm.pca <- train_caret.glm(mo.train, c("center", "scale", "nzv", "pca"))

  # results
  panel.stats.glm <- evaluate_regression_model(panel.test$y, predict(panel.fit.glm, panel.test))
  panel.stats.glm.pca <- evaluate_regression_model(panel.test$y, predict(panel.fit.glm.pca, panel.test))
  panel.stats.rf <- evaluate_regression_model(panel.test$y, predict(panel.fit.rf, panel.test))
  panel.stats.rf.pca <- evaluate_regression_model(panel.test$y, predict(panel.fit.rf.pca, panel.test))

  mo.stats.glm <- evaluate_regression_model(mo.test$y, predict(mo.fit.glm, mo.test))
  mo.stats.glm.pca <- evaluate_regression_model(mo.test$y, predict(mo.fit.glm.pca, mo.test))
  mo.stats.rf <- evaluate_regression_model(mo.test$y, predict(mo.fit.rf, mo.test))
  mo.stats.rf.pca <- evaluate_regression_model(mo.test$y, predict(mo.fit.rf.pca, mo.test))
  # output
  return(list(
    panel = list(
      "panel.fit" = panel.fit.glm, "panel.stats" = panel.stats.glm,
      "panel.fit.pca" = panel.fit.glm.pca, "panel.stats.pca" = panel.stats.glm.pca,
      "panel.test" = panel.test,
      "panel.fit.rf" = panel.fit.rf, "panel.stats.rf" = panel.stats.rf,
      "panel.fit.rf.pca" = panel.fit.rf.pca, "panel.stats.rf.pca" = panel.stats.rf.pca
    ),
    mo = list(
      "mo.fit" = mo.fit.glm, "mo.stats" = mo.stats.glm,
      "mo.fit.pca" = mo.fit.glm.pca, "mo.stats.pca" = mo.stats.glm.pca,
      "mo.test" = mo.test,
      "mo.fit.rf" = mo.fit.rf, "mo.stats.rf" = mo.stats.rf,
      "mo.fit.rf.pca" = mo.fit.rf.pca, "mo.stats.rf.pca" = mo.stats.rf.pca
    )
  ))
}

tic(msg = "Modelling", quiet = FALSE, func.tic = my.msg.tic)
# remove missing
dr <- dr[!is.na(dr$value)][!is.na(sample_id)][!is.na(column_name)]
dr$column_name <- gsub("/", "-", dr$column_name)
#candidates <- names(which(table(dr[column_name != "untreated"]$column_name) > as.numeric(args[2])))
candidates <- names(table(dr[column_name != "untreated"]$column_name))
candidates <- candidates[!(candidates %in% "")]

# assign a unique identifier to the modelling run
run.code <- paste0(sample(1:1e2, 1), sample(letters, 1))
outdir <- paste0(run.code, "_", dset, "_caretRes")
if (!dir.exists(file.path(p.out, outdir))) {
  dir.create(file.path(p.out, outdir))
}  

# start the parallelization
cl <- parallel::makeForkCluster(15)
doParallel::registerDoParallel(cl)
results <- foreach(drug = candidates) %dopar% {
  r <- run_caret(dat, dr, drugName = drug)
  saveRDS(r, file = file.path(p.out, outdir, paste0(drug, ".caret.RDS")))
  return(r)
}
parallel::stopCluster(cl)
toc(quiet = FALSE, func.toc = my.msg.toc)

## assemble model stats
names(results) <- candidates
dt <- do.call(
  rbind,
  lapply(
    names(results),
    function(drug) {
      x <- results[[drug]]
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
saveRDS(dt, file = file.path(p.out, outdir, "caret.stats.RDS"))

## Estimate feature importance glmnet
plan(multisession, workers = 4)
list.files(file.path(p.out, outdir),
           pattern = ".caret.RDS",
           full.names = T
) %>%
  future_map(
    .,
    ~ caret::varImp(readRDS(.x)$mo$mo.fit)$importance
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
message(date(), "\nFinished!")
