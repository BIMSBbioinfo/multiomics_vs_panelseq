args <- commandArgs(trailingOnly = TRUE)

# {SETUP}
## Paths
p.parent.dir <- getwd()
dset <- args[1] ## Either CCLE or PDX
p.data <- args[2] ## Path to a folder that contains "prepared" data
p.out <- args[3] ## Path to a folder where the models and their summary statistics will be written
## libraries & functions
if (!require("foreach")) install.packages("foreach") else library(foreach)
if (!require("data.table")) install.packages("data.table") else library(data.table)
source("src/F_auxiliary.R")
## read in data
dat <- readRDS(file.path(p.data, "prepared", paste0("data_", dset, ".RDS")))
dr <- data.table::fread(file.path(p.data, "Raw", dset, "drug_response.tsv")) # get drug response data

# {MAIN}
# df: samples on rows, features on columns, includes outcome variable 'y'
train_caret.rf <- function(df, ppOpts = c("center", "scale"), tgrid = NULL, folds = 5, reps = 1) {
  require(caret)
  set.seed(1234)
  model_caret <- train(y ~ .,
    data = df,
    method = "ranger",
    trControl = trainControl(
      method = "repeatedcv", number = folds,
      verboseIter = T, repeats = reps
    ),
    tuneGrid = tgrid,
    preProcess = ppOpts,
    num.trees = 500,
    importance = "permutation",
    num.threads = 2
  )
  return(model_caret)
}

train_caret.glm <- function(df, ppOpts = c("center", "scale"), tgrid = NULL, folds = 5, reps = 1) {
  require(caret)
  set.seed(1234)
  model_caret <- train(y ~ .,
    data = df,
    method = "glmnet",
    trControl = trainControl(
      method = "repeatedcv", number = folds,
      verboseIter = T, repeats = reps
    ),
    preProcess = ppOpts,
    tuneGrid = tgrid
  )
  return(model_caret)
}

# dat: prepared data including omics + drug response
# dr: drug response table
# drugName: name of drug to be analysed
# ppOpts: preprocessing options to be used during modeling  (e.g. scale/center/nzv/pca)
# algorithm: "RF" for random forests or "GLM" for elastic nets using GLMNET 
run_caret <- function(dat, dr, drugName, ppOpts, algorithm, tgrid = NULL, folds = 5, reps = 1) {
  set.seed(1234)
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

  message(date(), " => preparing datasets")
  # prepare data for panel (mut + cnv)
  panel.train <- data.frame(
    do.call(cbind, lapply(
      dat.raw[c("mut.panel", "cnv.panel")],
      function(x) t(x[, train_samples])
    )),check.names = F)
  panel.train$y <- y.train
  
  panel.test <- data.frame(do.call(cbind, lapply(
    dat.raw[c("mut.panel", "cnv.panel")],
    function(x) t(x[, test_samples])
  )), check.names = F)
  panel.test$y <- y.test
  
  # prepare data for mut+cnv+gex features (mo=multiomics) -------------------
  mo.train <- data.frame(do.call(cbind, lapply(
    dat.raw[c("mut.panel", "cnv.panel", "gex")],
    function(x) t(x[, train_samples])
  )), check.names = F)
  mo.train$y <- y.train
  
  mo.test <- data.frame(do.call(cbind, lapply(
    dat.raw[c("mut.panel", "cnv.panel", "gex")],
    function(x) t(x[, test_samples])
  )), check.names = F)
  mo.test$y <- y.test
  
  # build models for both panel and multiomics
  if(algorithm == 'RF') {
    message(date(), " => modeling for panel features")
    panel.fit <- train_caret.rf(df = panel.train, ppOpts = ppOpts, tgrid = tgrid, folds = folds, reps = reps)
    message(date(), " => modeling for multiomics features")
    mo.fit <- train_caret.rf(df = mo.train, ppOpts = ppOpts, tgrid = tgrid, folds = folds, reps = reps)
  } else if(algorithm == 'GLM') {
    message(date(), " => modeling for panel features")
    panel.fit <- train_caret.glm(df = panel.train, ppOpts = ppOpts, tgrid = tgrid, folds = folds, reps = reps)
    message(date(), " => modeling for multiomics features")
    mo.fit <- train_caret.glm(df = mo.train, ppOpts = ppOpts, tgrid = tgrid, folds = folds, reps = reps)
  }

  # extract stats
  panel.stats <- evaluate_regression_model(panel.test$y, stats::predict(panel.fit, panel.test))
  mo.stats <- evaluate_regression_model(mo.test$y, stats::predict(mo.fit, mo.test))
  
  stats <- rbind(panel.stats, mo.stats)
  stats$drugName <- drugName
  stats$type <- c("panel", "multiomics")
  stats$model <- algorithm
  stats$ppOpts <- paste(ppOpts, collapse = '+')
  stats$total_sample_count <- length(selected)
  stats$training_sample_count <- length(train_samples)
  stats$testing_sample_count <- length(test_samples)
  
  return(list('panel.fit' = panel.fit, 'mo.fit' = mo.fit, 
              'stats' = stats))
}

# remove missing
dr <- dr[!is.na(dr$value)][!is.na(sample_id)][!is.na(column_name)]
dr$column_name <- gsub("/", "-", dr$column_name)
# pick drugs with treated on at least 100 samples 
candidates <- dr[!column_name %in% c('untreated', ''),length(unique(sample_id)),by = column_name][V1 > 100]$column_name
#candidates <- c('TPCA-1', 'VENETOCLAX', 'METHOTREXATE', 'I-BET-762', 'LEFLUNOMIDE')

message(date()," => Modelling for ",length(candidates)," drugs")
outdir <- dset 
# create a subfolder on the target path for each dataset
if (!dir.exists(file.path(p.out, outdir))) {
  dir.create(file.path(p.out, outdir))
}  

# Define hyperparameter optimisation tuning grids for different algorithms
# random forests
tgrid_rf <- expand.grid(
   .mtry = seq(from = 10, to = 30, 10),
   .splitrule = "variance",
   .min.node.size = c(10, 20)
)

# using default tgrid for GLMNet
tgrid_glm <- expand.grid(
  .alpha = seq(0, 1, length = 10),
  .lambda = seq(0.0001, 1, length = 20)
)

cl <- parallel::makeForkCluster(30)
doParallel::registerDoParallel(cl)
foreach(drug = candidates) %dopar% {
  f <- file.path(p.out, outdir, paste0(drug, ".caret.RDS"))
  if(!file.exists(f)) {
    # results using random forests
    rf <- run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv"), 
                    algorithm = 'RF', tgrid = tgrid_rf, folds = 5, reps = 1)
    rf_pca <- run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv", "pca"), 
                        algorithm = 'RF', tgrid = tgrid_rf, folds = 5, reps = 2)
    
    glm <- run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv"), 
                     algorithm = 'GLM', tgrid = tgrid_glm, folds = 5, reps = 2)
    glm_pca <- run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv", "pca"), 
                         algorithm = 'GLM', tgrid = tgrid_glm, folds = 5, reps = 2)
    
    results <- list('rf' = rf, 'rf_pca' = rf_pca, 'glm' = glm, 'glm_pca' = glm_pca)
    saveRDS(results, file = f)
  }
}
parallel::stopCluster(cl)

message(date(), " => Finished!")
