#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)


# {SETUP}
## libraries & functions
library(foreach)
library(tictoc)
library(data.table)
utility_script <- "/data/local/buyar/arcas/uyar_et_al_multiomics_deeplearning/src/common_functions.R" # utility script
source(utility_script)
# time functions
my.msg.tic <- function(tic, msg) {
  if (is.null(msg) || is.na(msg) || length(msg) == 0)
  {
    outmsg <- paste0("Finished\nElapsed time: ", lubridate::seconds_to_period(round(toc - tic, 0)))
  }
  else
  {
    outmsg <- paste0("Starting ", msg, "...")
  }
}
my.msg.toc <- function(tic, toc, msg, info) {
  if (is.null(msg) || is.na(msg) || length(msg) == 0)
  {
    outmsg <- paste0("Finished\nElapsed time: ", lubridate::seconds_to_period(round(toc - tic, 0)))
  }
  else
  {
    outmsg <- paste0(msg, 
                     " finished\nElapsed time: ", 
                     lubridate::seconds_to_period(round(toc - tic, 0)),
                     "\n")
  }
}


## Paths
dset <-as.character(args[1])
p.data <- "/local/abarano/Projects/DrugResponse/data"
p.out <- "/local/abarano/Projects/DrugResponse/Results"
dat <- readRDS(file.path(p.data, paste0("data_", dset, ".RDS")))
dr <- data.table::fread(file.path(p.data, "Raw", dset, "drug_response.tsv")) # get drug response data 

# {MAIN}
# df: samples on rows, features on columns, includes outcome variable 'y'
train_caret.rf <- function(df, ppOpts = c("center", "scale")) {  
  require(caret)
  tgrid <- expand.grid(
    .mtry =  seq(from = 10, to = round(sqrt(ncol(df))), 10),
    .splitrule = 'variance',
    .min.node.size = c(5, 10, 15)
  )
  cl <- parallel::makePSOCKcluster(2)
  doParallel::registerDoParallel(cl)
  model_caret <- train(y  ~ ., 
                       data = df,
                       method = "ranger",
                       trControl = trainControl(method = "repeatedcv", number = 5, 
                                                verboseIter = T, repeats = 3),
                       tuneGrid = tgrid,
                       preProcess = ppOpts,
                       num.trees = 500,
                       importance = "permutation", 
                       num.threads = 2)
  parallel::stopCluster(cl)
  return(model_caret)
}
train_caret.glm <- function(df, ppOpts = c("center", "scale")) {  
  require(caret)
  tgrid <- expand.grid(.alpha = seq(0, 1, length = 10),
                       .lambda = seq(0.0001, 1, length = 20)
                       )
  model_caret <- train(y  ~ ., 
                       data = df,
                       method = "glmnet",
                       trControl = trainControl(method = "repeatedcv", number = 5, 
                                                verboseIter = T, repeats = 3),
                       tuneGrid = tgrid,
                       preProcess = ppOpts)
  return(model_caret)
}
run_caret <- function(dat, dr, drugName) {
  selected <- intersect(dr[column_name == drugName]$sample_id, colnames(dat$mut))
  colData <- dr[column_name == drugName][match(selected, sample_id)]

  # change feature names to avoid overlaps from different layers
  dat.raw <- sapply(simplify = F, names(dat), function(x) {
    m <- dat[[x]][,selected]
    rownames(m) <- paste0(rownames(m), '_', x)
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
  panel.train <- data.frame(do.call(cbind, 
                                    lapply(dat.raw[c('mut.panel', 'cnv.panel')], 
                                           function(x) t(x[,train_samples])
                                           )
                                    ), 
                            check.names = F)
  panel.train$y <- y.train
  panel.test <- data.frame(do.call(cbind, lapply(dat.raw[c('mut.panel', 'cnv.panel')], 
                                                  function(x) t(x[,test_samples])
                                                 )
                                   ), 
                            check.names = F)
  panel.test$y <- y.test
  # panel.fit.rf <- train_caret.rf(panel.train, c("center", "scale", "nzv"))
  # panel.fit.rf.pca <- train_caret.rf(panel.train, c("center", "scale", "nzv", "pca"))
  panel.fit.glm <- train_caret.glm(panel.train, c("center", "scale", "nzv"))
  panel.fit.glm.pca <- train_caret.glm(panel.train, c("center", "scale", "nzv", "pca"))
  
  
  
  
  # compute results for mut+cnv+gex features (mo=multiomics) -------------------
  message(date(), " => processing multiomics")
  mo.train <- data.frame(do.call(cbind, lapply(dat.raw[c('mut.panel', 'cnv.panel', 'gex')], 
                                               function(x) t(x[,train_samples]))), 
                         check.names = F)
  mo.train$y <- y.train
  mo.test <- data.frame(do.call(cbind, lapply(dat.raw[c('mut.panel', 'cnv.panel', 'gex')], 
                                              function(x) t(x[,test_samples]))), 
                        check.names = F)
  mo.test$y <- y.test
  # mo.fit.rf <- train_caret.rf(mo.train, c("center", "scale", "nzv"))
  # mo.fit.rf.pca <- train_caret.rf(mo.train, c("center", "scale", "nzv", "pca"))
  mo.fit.glm <- train_caret.glm(mo.train, c("center", "scale", "nzv"))
  mo.fit.glm.pca <- train_caret.glm(mo.train, c("center", "scale", "nzv", "pca"))

  
  # results
  panel.stats.glm <- evaluate_regression_model(panel.test$y, predict(panel.fit.glm, panel.test))
  panel.stats.glm.pca <- evaluate_regression_model(panel.test$y, predict(panel.fit.glm.pca, panel.test))
  mo.stats.glm <- evaluate_regression_model(mo.test$y, predict(mo.fit.glm, mo.test))
  mo.stats.glm.pca <- evaluate_regression_model(mo.test$y, predict(mo.fit.glm.pca, mo.test))
  
  # output
  return(list(panel = list('panel.fit' = panel.fit.glm, 'panel.stats' = panel.stats.glm,
                           'panel.fit.pca' = panel.fit.glm.pca, 'panel.stats.pca' = panel.stats.glm.pca,
                           'panel.test' = panel.test),
              mo = list('mo.fit' = mo.fit.glm, 'mo.stats' = mo.stats.glm,
                        'mo.fit.pca' = mo.fit.glm.pca, 'mo.stats.pca' = mo.stats.glm.pca,
                        'mo.test' = mo.test)
              )
        )
}

tic(msg = "Modelling", quiet = FALSE, func.tic = my.msg.tic)
# remove missing
dr <- dr[!is.na(dr$value)][!is.na(sample_id)][!is.na(column_name)]
dr$column_name <- gsub("/", "-", dr$column_name)
candidates <- names(which(table(dr[column_name != 'untreated']$column_name) > as.numeric(args[2])))

run.code <- paste0(sample(1:1e2, 1), sample(letters, 1))
outdir <- paste0(run.code, "_caretRes")
if (!dir.exists(file.path(p.out, outdir))) {
  dir.create(file.path(p.out, outdir))
}
# start the parallelization
cl <- parallel::makeForkCluster(12)
doParallel::registerDoParallel(cl)
results <- foreach(drug = candidates) %dopar% {
  r <- run_caret(dat, dr, drugName = drug)
  saveRDS(r, file = file.path(p.out, outdir, paste0(drug, ".caret.RDS")))
  return(r)
}
parallel::stopCluster(cl)
# cl <- parallel::makeCluster(12)
# parallel::clusterExport(cl = cl, varlist = c('utility_script', 'dat', 'dr',
#                                              'run_caret',
#                                              'train_caret.rf', 'train_caret.glm',
#                                              'p.out', 'outdir'))
# results <- pbapply::pblapply(cl = cl, 
#                              candidates, 
#                              function(d) {
#                                  source(utility_script)
#                                  require(data.table)
#                                  r <- run_caret(dat, dr, drugName = d)
#                                  saveRDS(r, file = file.path(p.out, outdir, paste0(d, ".caret.RDS")))
#                                  return(r)
#                                  }
#                             )
# parallel::stopCluster(cl)
toc(quiet = FALSE, func.toc = my.msg.toc)


names(results) <- candidates
dt <- do.call(rbind, 
              lapply(names(results), 
                     function(drug) {
                         x <- results[[drug]]
                         dt <- rbind(x$panel$panel.stats, 
                                     x$panel$panel.stats.pca,
                                     x$mo$mo.stats,
                                     x$mo$mo.stats.pca)
                         dt$type <- rep(c("panel", "mo"), each = 2)
                         dt$pp <- rep(c("scale+nzv", "scale+nzv+pca"), 2)
                         dt$model <- rep("glm", 4)
                         dt$drug <- rep(drug, 4)
                         return(dt)
                       }
                    )
              )
saveRDS(dt, file = file.path(p.out, outdir, "caret.stats.RDS"))
message(date(), "\nFinished!")


