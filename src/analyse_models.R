# analyse models built for drug response
## libraries & functions
if (!require("data.table")) install.packages("data.table") else library(data.table)
library(pbapply)
source("src/F_auxiliary.R")

args <- commandArgs(trailingOnly = TRUE)

folder <- args[1] # folder containing caret.RDS files (modeling results)

# start the parallelization
message(date()," => assembling model stats")

files <- dir(folder, ".caret.RDS$", full.names = T)

pbo = pbapply::pboptions(type="txt")

message(date()," => collating stats")
cl <- parallel::makeCluster(30)
parallel::clusterExport(cl = cl, varlist = c('files'))
stats <- do.call(rbind, pbapply::pblapply(cl = cl, files, function(f) {
  x <- readRDS(f)
  do.call(rbind, lapply(x, function(y) y$stats))
}))
parallel::stopCluster(cl)
saveRDS(stats, file = file.path(folder, "caret.stats.RDS"))


message(date()," => collating variable importance metrics")
# go through each multiomics model and extract variable importance metrics
# memory intensive job is to read the RDS files, so, we don't keep the contents in memory
cl <- parallel::makeCluster(30)
parallel::clusterExport(cl = cl, varlist = c('files'))
imp <- pbapply::pblapply(cl = cl, files, function(f) {
  require(data.table)
  x <- readRDS(f)
  drug <- gsub(".caret.RDS$", "", basename(f))
  imp <- caret::varImp(x[['rf']][['mo.fit']])[['importance']]
  dt.rf <- data.table('feature' = rownames(imp), 'val' = imp$Overall)
  colnames(dt.rf)[2] <- drug
  # parse glm
  imp <- caret::varImp(x[['glm']][['mo.fit']])[['importance']]
  dt.glm <- data.table('feature' = rownames(imp), 'val' = imp$Overall)
  colnames(dt.glm)[2] <- drug
  list('rf' = dt.rf, 'glm' = dt.glm)
})
parallel::stopCluster(cl)

rf_imp <- lapply(imp, function(x) x[['rf']])
glm_imp <- lapply(imp, function(x) x[['glm']])

rf_imp <- Reduce(function(x, y) merge(x, y, all=TRUE, by = 'feature'), rf_imp)
rf_imp[is.na(rf_imp)] <- 0 # Convert NA to zero
write.table(rf_imp, file = file.path(folder, "feature_importance.RF.tsv"), sep = '\t',
            quote = F, row.names = F)

glm_imp <- Reduce(function(x, y) merge(x, y, all=TRUE, by = 'feature'), glm_imp)
glm_imp[is.na(glm_imp)] <- 0 # Convert NA to zero
write.table(glm_imp, file = file.path(folder, "feature_importance.GLM.tsv"), sep = '\t',
            quote = F, row.names = F)

message(date(), " => finished assembling the model statistics")
