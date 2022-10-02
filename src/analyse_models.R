# analyse models built for drug response
## libraries & functions
if (!require("data.table")) install.packages("data.table") else library(data.table)
source("src/F_auxiliary.R")

args <- commandArgs(trailingOnly = TRUE)

folder <- args[1] # folder containing caret.RDS files (modeling results)

# start the parallelization
message(date()," => assembling model stats")

files <- dir(folder, ".caret.RDS$", full.names = T)

# read results files (includes fit models and statistics)
cl <- parallel::makeCluster(5)
parallel::clusterExport(cl = cl, varlist = c('files'))
results <- pbapply::pblapply(cl = cl, files, function(f) {
  x <- readRDS(f)
})
parallel::stopCluster(cl)
names(results) <- gsub(".caret.RDS$", "",basename(files))

# get stats
stats <- do.call(rbind, lapply(results, function(x) {
  do.call(rbind, lapply(x, function(y) y$stats))
}))
saveRDS(stats, file = file.path(folder, "caret.stats.RDS"))

# go through each multiomics model and extract variable importance metrics
get_varImp <- function(results, algorithm, type = 'mo.fit') {
  imp <- sapply(simplify = F, names(results), function(drug) {
    message(drug)
    imp <- caret::varImp(results[[drug]][[algorithm]][[type]])[['importance']]
    dt <- data.table('feature' = rownames(imp), 'val' = imp$Overall)
    colnames(dt)[2] <- drug
    return(dt)
  })
  imp <- Reduce(function(x, y) merge(x, y, all=TRUE, by = 'feature'), imp)
  # convert NA values to zero
  imp[is.na(imp)] <- 0
  return(imp)
}

# for RF:
rf_imp <- get_varImp(results, 'rf', 'mo.fit')
write.table(rf_imp, file = file.path(folder, "feature_importance.RF.tsv"), sep = '\t',
            quote = F, row.names = F)

# for GLM
glm_imp <- get_varImp(results, 'glm', 'mo.fit')
write.table(glm_imp, file = file.path(folder, "feature_importance.GLM.tsv"), sep = '\t',
            quote = F, row.names = F)

message(date(), " => finished assembling the model statistics")