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
    .mtry = seq(from = 5, to = round(sqrt(ncol(df))), 5),
    .splitrule = c("variance", "extratrees", "maxstat", "beta"),
    .min.node.size = seq(5, 25, 5)
  )
  cl <- parallel::makePSOCKcluster(4)
  doParallel::registerDoParallel(cl)
  set.seed(1234)
  model_caret <- train(y ~ .,
    data = df,
    method = "ranger",
    trControl = trainControl(
      method = "repeatedcv", number = 5,
      verboseIter = T, repeats = 5
    ),
    tuneGrid = tgrid,
    preProcess = c("center", "scale", "nzv"),
    num.trees = 1500,
    importance = "permutation",
    num.threads = 2
  )
  parallel::stopCluster(cl)
  return(model_caret)
}

train_caret.glm <- function(df) {
  require(caret)
  tgrid <- expand.grid(
    .alpha = seq(0, 1, length = 10),
    .lambda = seq(0.0001, 1, length = 20)
  )
  cl <- parallel::makePSOCKcluster(8)
  doParallel::registerDoParallel(cl)
  model_caret <- train(y ~ .,
    data = df,
    method = "glmnet",
    trControl = trainControl(
      method = "repeatedcv", number = 5,
      verboseIter = T, repeats = 3
    ),
    preProcess = c("center", "scale", "nzv"),
    tuneGrid = tgrid
  )
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

  mut <- t(dat$mut[, selected])
  gex <- t(dat$gex_gs[, selected])

  train_samples <- sample(selected, round(length(selected) * 0.7))
  test_samples <- setdiff(selected, train_samples)

  y.train <- drugs[variable == drugName][match(train_samples, name)]$value
  y.test <- drugs[variable == drugName][match(test_samples, name)]$value

  # compute results for mut features
  message(date(), " => processing panel")
  panel.train <- cbind(data.frame(mut[train_samples, ]), data.frame("y" = y.train))
  panel.test <- cbind(data.frame(mut[test_samples, ]), data.frame("y" = y.test))
  panel.fit.glm <- train_caret.glm(panel.train)
  panel.fit <- train_caret(panel.train)

  # compute results for mut+gex features (mo=multiomics)
  message(date(), " => processing multiomics")
  mo.train <- cbind(
    mut[train_samples, ], gex[train_samples, ],
    data.frame("y" = y.train)
  )
  mo.test <- cbind(
    mut[test_samples, ], gex[test_samples, ],
    data.frame("y" = y.test)
  )
  mo.fit.glm <- train_caret.glm(mo.train)
  mo.fit <- train_caret(mo.train)
  
  panel.stats <- evaluate_regression_model(panel.test$y, predict(panel.fit, panel.test))
  mo.stats <- evaluate_regression_model(mo.test$y, predict(mo.fit, mo.test))
  panel.stats.glm <- evaluate_regression_model(panel.test$y, predict(panel.fit.glm, panel.test))
  mo.stats.glm <- evaluate_regression_model(mo.test$y, predict(mo.fit.glm, mo.test))

  stats <- rbind(panel.stats, mo.stats, panel.stats.glm, mo.stats.glm)
  stats$type <- rep(c("panel", "multiomics"), 2)
  stats$model <- c(rep("rf", 2), rep("glm", 2))
  stats$total_sample_count <- length(selected)
  stats$training_sample_count <- length(train_samples)
  stats$testing_sample_count <- length(test_samples)

  # output
  rt <- list(
    panel = list(
      "panel.fit" = panel.fit, "panel.stats" = panel.stats,
      "panel.test" = panel.test,
      "panel.fit.glm" = panel.fit.glm, "panel.stats.glm" = panel.stats.glm
    ),
    mo = list(
      "mo.fit" = mo.fit, "mo.stats" = mo.stats,
      "mo.test" = mo.test,
      "mo.fit.glm" = mo.fit.glm, "mo.stats.glm" = mo.stats.glm
    ))
  saveRDS(rt, file = paste0("data/beatAML/",drugName, ".caret.RDS"))
  
  return(stats)
}

message(date(), "=> Started modelling")

drugs <- dat$drugs
candidates <- as.character(drugs[!is.na(value), length(name), by = variable][V1 > 100]$variable)

if (!dir.exists(file.path("data/beatAML"))) {
  dir.create(file.path("data/beatAML"))
}  

cl <- parallel::makeCluster(8)
parallel::clusterExport(cl = cl, varlist = c(
  "dat", "drugs", "evaluate_regression_model",
  "run_caret", "train_caret", "train_caret.glm"
))
results <- do.call(rbind, pbapply::pblapply(cl = cl, candidates, function(d) {
  require(data.table)
  r <- run_caret(dat, drugs, drugName = d)
  r$drug <- d
  return(r)
}))
parallel::stopCluster(cl)

# save stats, tables, figures
saveRDS(results, file = "beatAML.stats.RDS")

# write table
write.table(dcast(results, drug + total_sample_count + training_sample_count + testing_sample_count ~ type + model,
  value.var = c("RMSE", "COR", "Rsquare")
),
file = "beatAML.stats.tsv", sep = "\t", quote = F
)

# make summary figure
dt <- dcast.data.table(results, drug ~ type + model, value.var = "Rsquare")
dt$improvement_glm <- dt$multiomics_glm - dt$panel_glm
dt$improvement_rf <- dt$multiomics_rf - dt$panel_rf
for (md in c("rf", "glm")) {
  require(dplyr)
  da <- dt %>% dplyr::select(drug, ends_with(md))
  names(da) <- gsub(pattern = paste0("_", md), replacement = "", x = names(da))
  p1 <- ggplot(
    da,
    aes(x = panel, y = multiomics)
  ) +
    geom_point(aes(color = improvement), size = 3) +
    geom_abline(slope = 1) +
    coord_fixed() +
    lims(x = c(0, 0.5), y = c(0, 0.5)) +
    theme_bw(base_size = 12) +
    scale_color_gradient2(low = "black", mid = "gray", high = "red") +
    labs(color = "Multiomics\nimprovement") +
    theme(legend.position = "top")

  p2 <- results[results$model == md] %>%
    ggboxplot(x = "type", y = "Rsquare", add = "jitter") +
    stat_compare_means(
      paired = T, method.args = list("alternative" = "greater"),
      label.x = 2.2, label.y = 0.3
    ) +
    theme_bw(base_size = 12) +
    theme(axis.title.y = element_blank()) +
    coord_flip()

  p <- cowplot::plot_grid(p1, p2,
    ncol = 1, rel_heights = c(3, 1)
  )

  ggsave(filename = paste0("beatAML.", md, ".plot.pdf"), plot = p, width = 5, height = 7)
}

#Variable importance
## Estimate feature importance glmnet
plan(multisession, workers = 4)
for (v in c("mo.fit","mo.fit.glm")) {
  list.files(file.path("data/beatAML/"),
             pattern = ".caret.RDS",
             full.names = T
  ) %>%
    furrr::future_map(
      .,
      ~ caret::varImp(readRDS(.x)$mo[[v]])$importance
    ) -> l.varimp
  names(l.varimp) <- list.files(file.path("data/beatAML/"),
                                pattern = ".caret.RDS",
                                full.names = F
  ) %>%
    basename() %>%
    stringr::str_remove(., ".caret.RDS")
  purrr::map2(
    l.varimp,
    names(l.varimp),
    ~ .x %>%
      tibble::as_tibble(rownames = "Feature") %>%
      mutate(drug = .y)
  ) %>%
    purrr::reduce(., rbind) -> t.varimp
  write_tsv(t.varimp, file.path("data/beatAML/",v, "varImp.tsv"))
}

message(date(), "=> Finished!")
