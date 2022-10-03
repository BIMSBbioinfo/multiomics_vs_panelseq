# analyse models built for drug response
## libraries & functions
if (!require("data.table")) install.packages("data.table") else library(data.table)
library(pbapply)
source("src/F_auxiliary.R")

args <- commandArgs(trailingOnly = TRUE)

folder <- args[1] # folder containing caret.RDS files (modeling results)
nCores <- args[2] # number of cores to use for parallelisation
# start the parallelization
message(date()," => assembling model stats")

files <- dir(folder, ".caret.RDS$", full.names = T)

pbo = pbapply::pboptions(type="txt")

message(date()," => collating stats")
cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl = cl, varlist = c('files'))
stats <- do.call(rbind, pbapply::pblapply(cl = cl, files, function(f) {
  x <- readRDS(f)
  do.call(rbind, lapply(x, function(procedure) {
    do.call(rbind, lapply(procedure, function(run) {
      run$stats
    }))
  }))
}))
parallel::stopCluster(cl)

saveRDS(stats, file = file.path(folder, "caret.stats.RDS"))

message(date()," => collating variable importance metrics")
# go through each multiomics model and extract variable importance metrics
# memory intensive job is to read the RDS files, so, we don't keep the contents in memory

algorithms <- grep('pca', names(x), invert = T, value = T)
cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl = cl, varlist = c('files', 'algorithms'))
imp <- pbapply::pblapply(cl = cl, files, function(f) {
  require(data.table)
  x <- readRDS(f)
  drug <- gsub(".caret.RDS$", "", basename(f))
  # extract variable importance metrics for the first "run" 
  # for each method (excluding ppOpts with PCA)
  # algorithm => run => type (multiomics)
  L <- sapply(simplify = F, algorithms, function(algorithm) {
    imp <- caret::varImp(x[[algorithm]][[1]][['mo.fit']])[['importance']]
    dt <- data.table('feature' = rownames(imp), 'val' = imp$Overall)
    colnames(dt)[2] <- drug
    return(dt)
  })
})
parallel::stopCluster(cl)

imp_dt <- sapply(simplify = F, algorithms, function(a) {
  # get per-drug importance metrics for each algorithm
  l <- lapply(imp, function(x) x[[a]])
  # merge them into a single table
  M <- Reduce(function(x, y) merge(x, y, all=TRUE, by = 'feature'), l)
  M[is.na(M)] <- 0 # Convert NA to zero
  return(M)
})

# save tables to file
lapply(names(imp_dt), function(a) {
  write.table(imp_dt[[a]], file = file.path(folder, 
                                            paste0("feature_importance.", a, ".tsv")), sep = '\t',
              quote = F, row.names = F)
})

message(date(), " => finished assembling the model statistics")
