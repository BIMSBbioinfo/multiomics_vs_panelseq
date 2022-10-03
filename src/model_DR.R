args <- commandArgs(trailingOnly = TRUE)

# {SETUP}
## Paths
p.parent.dir <- getwd()
dset <- args[1] ## Either CCLE or PDX
p.data <- args[2] ## Path to a folder that contains "prepared" data
p.out <- args[3] ## Path to a folder where the models and their summary statistics will be written
repeatModeling <- as.numeric(args[4]) # Number of times the modeling procedure should be repeated with resampled train/test splits
nCores <- as.numeric(args[5]) #number cores to use for parallelisation

## libraries & functions
library(data.table)
library(caret)
library(pbapply)

utility_script <- "src/F_auxiliary.R"
source(utility_script)
## read in data
dat <- readRDS(file.path(p.data, "prepared", paste0("data_", dset, ".RDS")))
dr <- data.table::fread(file.path(p.data, "Raw", dset, "drug_response.tsv")) # get drug response data


train_caret <- function(df, algorithm, ppOpts = c("center", "scale"), tgrid = NULL, folds = 5, reps = 1, ...) {
  require(caret)
  set.seed(1234)
  model_caret <- train(y ~ .,
                       data = df,
                       method = algorithm,
                       trControl = trainControl(
                         method = "repeatedcv", number = folds,
                         verboseIter = T, repeats = reps
                       ), 
                       preProcess = ppOpts,
                       tuneGrid = tgrid,
                       ...
  )
  return(model_caret)
}

# prepare panel/multiomics datasets to be modeled 
# consolidate datasets, apply train/test splits 
# dat: omics data list
# dat: prepared data including omics + drug response
# dr: drug response table
# drugName: name of drug to be analysed
# setSeed: TRUE/FALSE whether to set a seed. 
prepareDataForModeling <- function(dat, dr, drugName, setSeed) {
  set.seed(setSeed)
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
  # prepare data for panel (mut.panel and/or cnv.panel)
  panel_layers <- grep('panel', names(dat.raw), value = T)
  panel.train <- data.frame(
    do.call(cbind, lapply(
      dat.raw[panel_layers],
      function(x) t(x[, train_samples])
    )),check.names = F)
  panel.train$y <- y.train
  
  panel.test <- data.frame(do.call(cbind, lapply(
    dat.raw[panel_layers],
    function(x) t(x[, test_samples])
  )), check.names = F)
  panel.test$y <- y.test
  
  # prepare data for panel features + gex (mo=multiomics) -------------------
  mo.train <- data.frame(do.call(cbind, lapply(
    dat.raw[c(panel_layers, "gex")],
    function(x) t(x[, train_samples])
  )), check.names = F)
  mo.train$y <- y.train
  
  mo.test <- data.frame(do.call(cbind, lapply(
    dat.raw[c(panel_layers, "gex")],
    function(x) t(x[, test_samples])
  )), check.names = F)
  mo.test$y <- y.test
  return(list('panel.train' = panel.train, 'panel.test' = panel.test, 
              'mo.train' = mo.train, 'mo.test' = mo.test, 
              'total_sample_count' = length(selected), 
              'training_sample_count' = length(train_samples), 
              'testing_sample_count' = length(test_samples)))
}
 
# dat: prepared data including omics + drug response
# dr: drug response table
# drugName: name of drug to be analysed
# ppOpts: preprocessing options to be used during modeling  (e.g. scale/center/nzv/pca)
# algorithm: "RF" for random forests or "GLM" for elastic nets using GLMNET 
# tgrid: tuning grid to use for the corresponding algorithm
# folds: number of folds for cross-validation
# reps: number of repeates for repeating the cross-validation procedure
# setSeed: seed value to set while sampling train/test samples
# runId: a numerical value that represents the current run
# ... : additional options to pass to train_caret -> train
run_caret <- function(dat, dr, drugName, ppOpts, algorithm, tgrid = NULL, folds = 5, reps = 1, runId, setSeed, ...) {
  
  ds <- prepareDataForModeling(dat, dr, drugName, setSeed = setSeed)
  panel.train <- ds$panel.train
  panel.test <- ds$panel.test
  mo.train <- ds$mo.train
  mo.test <- ds$mo.test
  
  # build models for both panel and multiomics
  message(date(), " => modeling for panel features")
  panel.fit <- train_caret(df = panel.train, algorithm = algorithm, ppOpts = ppOpts, 
                           tgrid = tgrid, folds = folds, reps = reps, ...)
  message(date(), " => modeling for multiomics features")
  mo.fit <- train_caret(df = mo.train, algorithm = algorithm, ppOpts = ppOpts, tgrid = tgrid, 
                           folds = folds, reps = reps, ...)
  # extract stats
  panel.stats <- evaluate_regression_model(panel.test$y, stats::predict(panel.fit, panel.test))
  mo.stats <- evaluate_regression_model(mo.test$y, stats::predict(mo.fit, mo.test))
  
  stats <- rbind(panel.stats, mo.stats)
  stats$drugName <- drugName
  stats$type <- c("panel", "multiomics")
  stats$model <- algorithm
  stats$ppOpts <- paste(ppOpts, collapse = '+')
  stats$run <- runId
  stats$seed <- setSeed
  stats$total_sample_count <- ds$total_sample_count
  stats$training_sample_count <- ds$training_sample_count
  stats$testing_sample_count <- ds$testing_sample_count
  
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

#in sample train/test splits, set seed only if the modeling is not repeated with different train/test splits
# generate seeds to be used in modeling runs
# the seed is always set to 1234 for the first run
seeds <- c(1234, sample(1:2^10, repeatModeling-1))  
cv_reps <- 1 #number of repetitions for cross-validation

pbo = pbapply::pboptions(type="txt")

cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl, varlist = c('candidates', 'seeds', 'cv_reps', 
                                        'dat', 'dr', 'repeatModeling', 
                                        'p.out', 'outdir',
                                        'run_caret', 'prepareDataForModeling', 
                                        'train_caret', 'utility_script'))
pbapply::pblapply(cl = cl, candidates, function(drug) { 
  source(utility_script)
  require(data.table)
  f <- file.path(p.out, outdir, paste0(drug, ".caret.RDS"))
  if(!file.exists(f)) {
    # results using random forests
    rf <- lapply(1:repeatModeling, function(i) {
      run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv"),
                algorithm = 'ranger', tgrid = NULL, folds = 5, reps = cv_reps, runId = i, 
                setSeed = seeds[i], num.threads = 2, importance = 'permutation')
    })
    rf_pca <-  lapply(1:repeatModeling, function(i) {
      run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv", "pca"), 
                algorithm = 'ranger', tgrid = NULL, folds = 5, reps = cv_reps, runId = i, 
                setSeed = seeds[i], num.threads = 2, importance = 'permutation')
    })
    glm <- lapply(1:repeatModeling, function(i) {
      run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv"), 
                algorithm = 'glmnet', tgrid = NULL, folds = 5, reps = cv_reps, runId = i, setSeed = seeds[i])
    })
    glm_pca <-  lapply(1:repeatModeling, function(i) {
      run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv", "pca"), 
                algorithm = 'glmnet', tgrid = NULL, folds = 5, reps = cv_reps, runId = i, setSeed = seeds[i])
    })
    
    svm <- lapply(1:repeatModeling, function(i) {
      run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv"), 
                algorithm = 'svmRadial', tgrid = NULL, folds = 5, reps = cv_reps, runId = i, setSeed = seeds[i])
    })
    svm_pca <- lapply(1:repeatModeling, function(i) {
      run_caret(dat = dat, dr = dr, drugName = drug, ppOpts = c("center", "scale", "nzv", "pca"), 
                algorithm = 'svmRadial', tgrid = NULL, folds = 5, reps = cv_reps, runId = i, setSeed = seeds[i])
    })
    
    results <- list('rf' = rf, 'rf_pca' = rf_pca, 'glm' = glm, 'glm_pca' = glm_pca, 
                    'svm' = svm, 'svm_pca' = svm_pca)
    saveRDS(results, file = f)
  }
})
parallel::stopCluster(cl)

message(date(), " => Finished!")
