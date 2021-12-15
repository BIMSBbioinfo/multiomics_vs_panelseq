#remove everything else
rm(list=ls()); gc()

#Load the libraries
library(caret)
library(tidyverse)
library(dplyr)
library(doParallel)
library(caretEnsemble)

#set the memory
memory.size(9.9e+10)
memory.limit(9.9e+10)

cores = 6 #number of the cores to be used

#start the parallelization
Mycluster <- makeCluster(cores) ; registerDoParallel(Mycluster)

# data specifics
l.d <- c("PDX_all")
d.type <- c("multiomics", "panelseq", "exomeseq")
set.seed(42)
run.code <- "70o"

p.d <- file.path("data/Results/")

#read the growth rate data
growth_rate <- read.delim("data/growth_rate.tsv") %>%
  dplyr::rename(sample_id = Model)%>%
  dplyr::select(sample_id,TimeToDouble)

# read the data
message("Reading in data...")
data.pca <- readRDS(file.path("data/Results", paste0("PCAembeddings_", run.code, ".RDS")))

# -------------------------------------------------------------------------------------------------
# prepare training&testing sets
message("Preparing training and testing partitions...")
training <- list()
testing <- list()
for(source in l.d) {
  for(type in d.type) {
    set.seed(312) #for reproducibility 
      tmp <- dplyr::inner_join(growth_rate, data.pca[[source]][[type]], by = "sample_id") %>% 
        as.data.frame() %>% na.omit()
      
      # get indices for 70% of the data set
      intrain <- createDataPartition(y = tmp[,2], p= 0.8)[[1]]
      
      # seperate test and training sets
      training[[source]][[type]] <- tmp[intrain,]
      testing[[source]][[type]] <- tmp[-intrain,]
  }
}
saveRDS(training,"data/growth_rate_Results/training_sets.RDS")
saveRDS(testing,"data/growth_rate_Results/testing_sets.RDS")
# -------------------------------------------------------------------------------------------------
tune_list <- list(
  svmLinear = caretModelSpec(method = "svmLinear2",
                             tuneGrid = expand.grid(cost = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25)),
                             importance = TRUE),
  gbm = caretModelSpec(method = "gbm",
                       tuneGrid = expand.grid(n.trees=c(100,200,300,400,500),
                                              interaction.depth=c(1,2),
                                              shrinkage=c(0.01,0.1),
                                              n.minobsinnode = 10)),
  ranger = caretModelSpec(method = "ranger",
                          tuneGrid = expand.grid(mtry = seq(2, 10, by = 1),
                                                 min.node.size = c(3, 5, 10),
                                                 splitrule = c("variance")),
                          importance = "impurity", 
                          num.threads = 1),
  glmnet = caretModelSpec(method = "glmnet",
                          tuneGrid = expand.grid(alpha = seq(0, 1, length = 30),
                                                 lambda = seq(0.0001, 1, length = 100)))
)
# -------------------------------------------------------------------------------------------------
models <- list()
for (source in l.d) {
  for(type in d.type){
    caretList(x = training[[source]][[type]][,-c(1,2)], 
              y = training[[source]][[type]]$TimeToDouble,
              trControl=trainControl(method="cv", number = 5,
                                     allowParallel = TRUE,
                                     savePredictions = "final",
                                     index = createResample(training[[source]][[type]]$sample_id),
                                     returnResamp = "final",
                                     search = "grid",
                                     verbose = TRUE),
              methodList = c("svmLinear2","rf","glmnet","gbm"),
              continue_on_fail = FALSE,
              tuneList=tune_list,
              verbose = TRUE) -> model.out
    models[[source]][[type]] <- model.out
  }
}

saveRDS(models, file.path("data/growth_rate_results/",  
                          paste0("growth_rate_Results_5kcv_", run.code, ".RDS")))
parallel::stopCluster(Mycluster)

# -------------------------------------------------------------------------------------------------
model_types <- c("svmLinear","ranger","glmnet","gbm")
summary_models <- function(model_list,validation_type = NA, test_list = NA,
                           model_types = NA){
  t.res <- dplyr::tibble()
  for (source in names(model_list)) {
    for (dt in names(model_list[[source]])) {
      #remove duplicates if there is any
      for (type in model_types) {
        model_list[[dt]][[type]] <- model_list[[source]][[dt]][[type]]
        
        #Extract the best tune
        tune <- model_list[[dt]][[type]]$bestTune
        tmp <- model_list[[dt]][[type]]$results
        for (i in names(tune)) {
          tmp <- tmp[tmp[[i]] == tune[[i]], ]
        }
        tmp <- tmp %>% 
          dplyr::select(RMSE, Rsquared, MAE) %>% 
          mutate(source = source,
                 dataset = dt,
                 model = type, 
                 validation = validation_type) %>% 
          dplyr::select(source, dataset, model, validation, everything())
        
        ## assemble
        t.res <- rbind(t.res, tmp) 
        if (!is.na(validation_type) && !is.na(test_list)) {
          ## test set validation
          tmp.2 <- tibble(obs = test_list[[source]][[dt]]$TimeToDouble, 
                          pred = model_list[[source]][[dt]][[type]] %>% 
                            predict(test_list[[source]][[dt]][, -c(1, 2)]))
          tmp.2_stats <- caret::postResample(tmp.2$pred, tmp.2$obs)
          tmp.2 <- rbind(tibble(), unname(tmp.2_stats))
          colnames(tmp.2) <- names(tmp.2_stats)
          tmp.2 <- tmp.2 %>% 
            dplyr::select(RMSE, Rsquared, MAE) %>% 
            mutate(source = source,
                   dataset = dt,
                   model = type, 
                   validation = "test partition") %>% 
            dplyr::select(source, dataset,  model, validation, RMSE, Rsquared, MAE)
          
          ## append
          t.res <- rbind(t.res, tmp.2)
          
        }
      }
    }
  }
  return(t.res)
}

# -------------------------------------------------------------------------------------------------
#summary results
testing_growth <- summary_models(models, "cv", test_list = testing,
                                 model_types = model_types)

#save all predictions
saveRDS(testing_growth,
        file.path("data/growth_rate_Results/growth_rate_predictions_all.RDS"))


# Calculates mean, sd,and se
tgrowth_sum <- testing_growth  %>%
  group_by(dataset,validation) %>%
  summarise( 
    n=n(),
    mean=mean(Rsquared),
    sd=sd(Rsquared)
  ) %>%
  mutate( se=sd/sqrt(n))

#plot the summary results
growth_plot <- tgrowth_sum %>%
  ggplot(aes(x= dataset,y=mean, fill=validation)) + 
  geom_bar(position = 'dodge', aes(x=dataset, y=mean), stat="identity") +
  geom_errorbar(aes(x=dataset, ymin=mean-se, ymax=mean+se), width=0.5,
                colour="black", size=1, position = position_dodge(0.9))+
  theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")+
  ggtitle("Growth Rate Predictions")+
  ylab("Mean of R-squared") + xlab("dataset")

#save the summary plot
pdf("misc/summarized_growth_Rate_prediction_plot.pdf", width = 20, height = 10)
growth_plot
dev.off()

#plot the models seperately
testing_growth %>%
  filter(Rsquared < 1) %>% 
  group_by(dataset, model) %>% 
  summarize(mean_rsqrd = mean(Rsquared),
            sd_rsqrd = sd(Rsquared)) 


  ggplot(data = , aes(x = dataset, y = mean_rsqrd)) + 
  geom_col(aes(fill = validation), position = "dodge", color = "grey30", width = 0.7) +
  facet_wrap(~ `model`) +
  geom_errorbar(aes(x=dataset, ymin=mean_rsqrd - sd_rsqrd, ymax=mean_rsqrd + sd_rsqrd,), width=0.5,
                colour="black", size=1, position = position_dodge(0.9))+
  geom_hline(yintercept = 0) +
  lims(y = c(-0.1, 0.5)) +
  labs(x = "Dataset", y = "Mean of Rsquared") +
   
  theme_bw(base_size = 14) + theme(aspect.ratio = 1.7, axis.text.x = element_text(angle = 45, hjust = 1)) 
