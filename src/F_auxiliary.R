## Auxiliary functions

# convenience
evaluate_regression_model <- function(y, y_hat) {
  require(caret)  
# Model performance metrics
  data.table::data.table(
    RMSE = RMSE(y_hat, y),
    Rsquare = R2(y_hat, y),
    COR = cor(y_hat, y)
  )
}
score_gene_set <- function(rankData, genes_up, gene_down = NULL) {
  scores <- singscore::simpleScore(rankData, upSet = genes_up)
  return(scores$TotalScore)
}

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


## Plotting funcitons
plot_pcaImp <- function(data_melted, base.size = 12) {
  (
    (ggplot(data = data_melted, aes(y = `scale+nzv+pca`, x = "")) + 
       geom_boxplot() + 
       lims(y = c(0, 0.5)) +
       theme_bw(base_size = base.size) + labs(x = "")
    ) +
      (ggplot(data = data_melted, aes(x = `scale+nzv`, y = `scale+nzv+pca`)) + 
         geom_point(aes(color = `scale+nzv+pca` - `scale+nzv`)) + 
         geom_abline(intercept = 0, slope = 1) + 
         lims(x = c(0, 0.5), y = c(0, 0.5)) +
         scale_color_gradientn(colours = terrain.colors(100)) +
         theme_bw(base_size = base.size) + labs(x = "", y = "", color = "PCA\nimprovement")
      ) +
      (plot_spacer()
      ) + 
      (ggplot(data = data_melted, aes(x = `scale+nzv`, y = "")) + 
         geom_boxplot() + 
         lims(x = c(0, 0.5)) +
         theme_bw(base_size = base.size) + labs(y = "")
      )
  ) + plot_layout(widths = c(0.2, 1),
                  heights = c(1, 0.2))
}
plot_moImp <- function(data_melted, base.size = 12) {
  (
    (ggplot(data = data_melted, aes(y = mo, x = "")) + 
       geom_boxplot() + 
       lims(y = c(0, 0.5)) +
       theme_bw(base_size = base.size) + labs(x = "")
    ) +
      (ggplot(data = data_melted, aes(x = panel, y = mo)) + 
         geom_point(aes(color = mo - panel)) + 
         geom_abline(intercept = 0, slope = 1) + 
         lims(x = c(0, 0.5), y = c(0, 0.5)) +
         scale_color_gradientn(colours = terrain.colors(100)) +
         theme_bw(base_size = base.size) + labs(x = "", y = "", color = "Multiomics\nimprovement")
      ) +
      (plot_spacer()
      ) + 
      (ggplot(data = data_melted, aes(x = panel, y = "")) + 
         geom_boxplot() + 
         lims(x = c(0, 0.5)) +
         theme_bw(base_size = base.size) + labs(y = "")
      )
  ) + plot_layout(widths = c(0.2, 1),
                  heights = c(1, 0.2))
}