# build drug response models for beatAML study
# test caret using ranger without dimension reduction

library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(caret)
library(ranger)
library(MetBrewer)
ggplot2::theme_set(ggpubr::theme_pubclean())

args <- commandArgs(trailingOnly = T)

beatAML_prepared_data <- args[1] 
 
dat <- readRDS(beatAML_prepared_data)

# define a function to train a random forest model using 
# a repeated-cross-validation procedure using caret. 
# df: samples on rows, features on columns, includes outcome variable 'y'
train_caret <- function(df) {
  require(caret)
  tgrid <- expand.grid(
    .mtry =  seq(from = 10, to = round(sqrt(ncol(df))), 10),
    .splitrule = 'variance',
    .min.node.size = c(5,10,15)
  )
  cl <- parallel::makePSOCKcluster(5)
  doParallel::registerDoParallel(cl)
  set.seed(1234)
  model_caret <- train(y  ~ ., 
                       data = df,
                       method = "ranger",
                       trControl = trainControl(method="repeatedcv", number = 5, 
                                                verboseIter = T, repeats = 5),
                       tuneGrid = tgrid,
                       preProcess = c('center', 'scale', 'nzv'), 
                       num.trees = 500,
                       importance = "permutation", 
                       num.threads = 2)
  parallel::stopCluster(cl)
  return(model_caret)
}

# function to evaluate a regression model's performance
evaluate_regression_model <- function(y, y_hat) {
  require(caret)  
  # Model performance metrics
  data.table::data.table(
    RMSE = RMSE(y_hat, y),
    Rsquare = R2(y_hat, y),
    COR = cor(y_hat, y)
  )
}

# This function builds two sets of models using a train/test split (70/30)
# The models will be trained on train split, stats will be generated on 
# the test split. 
# 1. using only mutation features (mut)
# 2. using mutation features (mut) + rna-seq-based gene-set scores (mut + gex)
# output: evaluation results of each model for each drug tested 
run_caret <- function(dat, drugs, drugName) {
  set.seed(1234)
  samples <- drugs[variable == drugName][!is.na(value)]$name
  selected <- Reduce(intersect, list(samples, colnames(dat$mut), colnames(dat$gex_gs)))
  
  mut <- t(dat$mut[,selected])
  gex <- t(dat$gex_gs[,selected])
  
  train_samples <- sample(selected, round(length(selected) * 0.7))
  test_samples <- setdiff(selected, train_samples)
  
  y.train <- drugs[variable == drugName][match(train_samples, name)]$value
  y.test <- drugs[variable == drugName][match(test_samples, name)]$value
  
  # compute results for mut features
  message(date(), " => processing panel")
  panel.train <- cbind(data.frame(mut[train_samples,]), data.frame('y' = y.train))
  panel.test <- cbind(data.frame(mut[test_samples,]), data.frame('y' = y.test))
  panel.fit <- train_caret(panel.train)

  # compute results for mut+gex features (mo=multiomics)
  message(date(), " => processing multiomics")
  mo.train <- cbind(mut[train_samples,], gex[train_samples,], 
                    data.frame('y' = y.train))
  mo.test <- cbind(mut[test_samples,], gex[test_samples,], 
                   data.frame('y' = y.test))
  mo.fit <- train_caret(mo.train)
  
  panel.stats <- evaluate_regression_model(panel.test$y, predict(panel.fit, panel.test))
  mo.stats <- evaluate_regression_model(mo.test$y, predict(mo.fit, mo.test))

  stats <- rbind(panel.stats, mo.stats)
  stats$type <- c('panel', 'multiomics')
  
  return(stats)
}

message(date(), "=> Started modelling")

drugs <- dat$drugs
candidates <- as.character(drugs[!is.na(value),length(name),by=variable][V1 > 100]$variable)

cl <- parallel::makeCluster(10)
parallel::clusterExport(cl = cl, varlist = c('dat', 'drugs', 'evaluate_regression_model',
                                             'run_caret', 'train_caret'))
results <- do.call(rbind, pbapply::pblapply(cl = cl, candidates, function(d) {
  require(data.table)
  r <- run_caret(dat, drugs, drugName = d)
  r$drug <- d
  return(r)
}))
parallel::stopCluster(cl)

# save stats, tables, figures
saveRDS(results, file = 'beatAML.stats.RDS')

# write table 
write.table(dcast(results, drug ~ type, value.var = c('RMSE', 'COR', 'Rsquare')), 
            file = 'beatAML.stats.tsv', sep = '\t', quote = F)

# make summary figure 

dt <- dcast.data.table(results, drug ~ type, value.var = 'Rsquare')
dt$improvement <- dt$multiomics - dt$panel
p1 <- ggplot(dt,
       aes(x = panel, y = multiomics)) +
  geom_point(aes(color = improvement), size = 3) + geom_abline(slope = 1) + 
  coord_fixed() + 
  lims(x = c(0, 0.5), y = c(0, 0.5)) + 
  theme_bw(base_size = 12) +
  scale_color_gradient2(low = 'black', mid = 'gray', high = 'red') +
  labs(color = "Multiomics\nimprovement")

p2 <- ggboxplot(results, x = 'type', y = 'Rsquare', add = 'jitter') + 
  stat_compare_means(paired = T, method.args = list('alternative' = 'greater')) +
  theme_bw(base_size = 12) 

p <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 1))

ggsave(filename = 'beatAML.plot.pdf', plot = p, width = 12, height = 6) #, width = 6, height = 6)
message(date(), "=> Finished!")




