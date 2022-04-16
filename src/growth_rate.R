#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)


# {SETUP}
## libraries & functions
library(foreach)
library(tictoc)
library(tidyr)
library(data.table)
utility_script <- "common_functions.R" # utility script
source(utility_script)


set.seed(123)
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
dset <- "PDX"
p.data <- "data/"
p.out <- "data/Results/"
dat <- readRDS(file.path(p.data, paste0("data_", dset, ".RDS")))
dr <- data.table::fread(file.path(p.data,"growth_rate.tsv"))  %>%
  dplyr::rename(sample_id = Model) %>%
  tidyr::drop_na() %>%
  dplyr::select(sample_id,TimeToDouble)
dr$TimeToDouble <- round(dr$TimeToDouble, digits=1)

dr <- dr %>%
  dplyr::mutate(duplicated = duplicated(TimeToDouble)) %>%
  dplyr::filter(!duplicated == TRUE)%>%
  dplyr::select(sample_id,TimeToDouble)

# {MAIN}
# df: samples on rows, features on columns, includes outcome variable 'y'
train_caret.rf <- function(df, ppOpts = c("center", "scale")) {  
  require(caret)
  tgrid <- expand.grid(
    .mtry =  seq(from = 10, to = round(sqrt(ncol(df))), 10),
    .splitrule = 'variance',
    .min.node.size = c(5, 10, 15)
  )
  cl <- parallel::makePSOCKcluster(20)
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
run_caret <- function(dat, dr) {
  selected <- intersect(dr$sample_id, colnames(dat$mut))
  colData <- dr[match(selected, sample_id)]
  
  # change feature names to avoid overlaps from different layers
  dat.raw <- sapply(simplify = F, names(dat), function(x) {
    m <- dat[[x]][,selected]
    rownames(m) <- paste0(rownames(m), '_', x)
    return(m)
  })
  train_samples <- sample(selected, round(length(selected) * 0.7))
  test_samples <- setdiff(selected, train_samples)
  
  y.train <- colData[match(train_samples, sample_id)]$TimeToDouble
  y.test <- colData[match(test_samples, sample_id)]$TimeToDouble
  
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
  panel.fit.rf <- train_caret.rf(panel.train, c("center", "scale","nzv"))
  panel.fit.rf.pca <- train_caret.rf(panel.train, c("center", "scale","nzv", "pca"))

  
  
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
  mo.fit.rf <- train_caret.rf(mo.train, c("nzv"))
  mo.fit.rf.pca <- train_caret.rf(mo.train, c( "nzv", "pca"))

  
  
  # results
  panel.stats.rf <- evaluate_regression_model(panel.test$y, predict(panel.fit.rf, panel.test))
  panel.stats.rf.pca <- evaluate_regression_model(panel.test$y, predict(panel.fit.rf.pca, panel.test))
  mo.stats.rf <- evaluate_regression_model(mo.test$y, predict(mo.fit.rf, mo.test))
  mo.stats.rf.pca <- evaluate_regression_model(mo.test$y, predict(mo.fit.rf.pca, mo.test))
  
  # output
  return(list(panelseq = list('panel.fit' = panel.fit.rf, 'panel.stats' = panel.stats.rf,
                           'panel.fit.pca' = panel.fit.rf.pca, 'panel.stats.pca' = panel.stats.rf.pca,
                           'panel.test' = panel.test),
              multiomics = list('mo.fit' = mo.fit.rf, 'mo.stats' = mo.stats.rf,
                        'mo.fit.pca' = mo.fit.rf.pca, 'mo.stats.pca' = mo.stats.rf.pca,
                        'mo.test' = mo.test)
  )
  )
}

tic(msg = "Modelling", quiet = FALSE, func.tic = my.msg.tic)
# remove missing
dr <- dr[!is.na(dr$TimeToDouble)][!is.na(sample_id)]

  # start the parallelization
  cl <- parallel::makeForkCluster(12)
  doParallel::registerDoParallel(cl)
  
  run.code <- paste0(sample(1:1e2, 1), sample(letters, 1))
  outdir <- paste0(run.code, "_caretRes")
  if (!dir.exists(file.path(p.out, outdir))) {
    dir.create(file.path(p.out, outdir))
  }
  results <-  run_caret(dat, dr)
  saveRDS(results, file = file.path(p.out, outdir, paste0(".caret.RDS")))
  
  parallel::stopCluster(cl)
  toc(quiet = FALSE, func.toc = my.msg.toc)
  
  dt <- do.call(rbind, 
                lapply(names(results)[1], 
                       function(i) {
                         x <- results
                         dt <- rbind(x$panelseq$panel.stats, 
                                     x$panelseq$panel.stats.pca,
                                     x$multiomics$mo.stats,
                                     x$multiomics$mo.stats.pca)
                         dt$type <- rep(c("panelseq", "multiomics"), each = 2)
                         dt$pp <- rep(c("scale+nzv", "scale+nzv+pca"), 2)
                         dt$model <- rep("rf", 4)
                         return(dt)
                       }
                )
  )



saveRDS(dt, file = file.path(p.out, outdir, "caret.stats.RDS"))
message(date(), "\nFinished!")

#plot the results
plot <- dt%>%
  ggplot(aes(x= type,y=Rsquare, fill=type))+
  geom_bar(position = 'dodge', aes(x=pp, y=Rsquare), stat="identity")+
  theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")+
  ylab("R-squared") + xlab("dataset")

#save the plot as pdf
pdf("data/growth_rate_results.pdf", width = 20, height = 10)
plot
dev.off()
