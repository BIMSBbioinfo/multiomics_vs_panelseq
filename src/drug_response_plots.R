
# ---------------------------helper function-------------------------------------------------
summary_models <-  function(model_list,validation_type = NA, test_list = NA,
                            model_types = NA, data_source) {
  library(tidyverse)
  t.res <- dplyr::tibble()
  for (source in names(model_list[[data_source]])) {
    for (drug in names(model_list[[data_source]][[source]])) {
      for (type in model_types) {
        #remove duplicates if there is any
        model_list[[data_source]][[source]][[drug]][[type]] <- model_list[[data_source]][[source]][[drug]][[type]]
        
        #Extract the best tune
        tune <- model_list[[data_source]][[source]][[drug]][[type]]$bestTune
        tmp <- model_list[[data_source]][[source]][[drug]][[type]]$results
        
        for (i in names(tune)) {
          tmp <- tmp[tmp[[i]] == tune[[i]], ]
        }
        tmp <- tmp %>% 
          dplyr::select(RMSE, Rsquared, MAE) %>% 
          mutate(source = source,
                 drug = drug,
                 model = type, 
                 validation = validation_type) %>% 
          dplyr::select(source, drug, model, validation, everything())
        
        ## assemble
        t.res <- rbind(t.res, tmp) 
        if (!is.na(validation_type) && !is.na(test_list)) {
          ## test set validation
          tmp.2 <- tibble(obs = test_list[[data_source]][[source]][[drug]]$value, 
                          pred = (model_list[[data_source]][[source]][[drug]][[type]] %>% 
                                    predict(test_list[[data_source]][[source]][[drug]][, -c(1, 2)])))
          tmp.2_stats <- caret::postResample(tmp.2$pred, tmp.2$obs)
          tmp.2 <- rbind(tibble(), unname(tmp.2_stats))
          colnames(tmp.2) <- names(tmp.2_stats)
          tmp.2 <- tmp.2 %>% 
            dplyr::select(RMSE, Rsquared, MAE) %>% 
            mutate(source = source,
                   drug = drug,
                   model = type, 
                   validation = "test partition") %>% 
            dplyr::select(source, drug,  model, validation, RMSE, Rsquared, MAE)
          
          ## append
          t.res <- rbind(t.res, tmp.2)
        }
      }
    }
  }
  # Calculates mean, sd,and se
  table <- t.res  %>%
    na.omit() %>%
    group_by(validation,source) %>%
    summarise( 
      n=n(),
      mean=mean(Rsquared),
      sd=sd(Rsquared)
    ) %>%
    mutate( se=sd/sqrt(n))
  
  #plot the results
  plot <- table %>%
    ggplot(aes(x= source,y=mean, fill=source)) + 
    geom_bar(position = 'dodge', aes(x=source, y=mean), stat="identity") +
    geom_errorbar(aes(x=source, ymin=mean-se, ymax=mean+se), width=0.5,
                  colour="black", size=1, position = position_dodge(0.9))+
    theme_bw()+
    theme(legend.position = "top",plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette="Set1")+
    ylim(c(0,0.5))+
    facet_wrap(~validation)+
    ylab("Mean of R-squared") + xlab("dataset")
  return(list(Summary = table, Plot = plot)) 
}

# -------------------------------------------------------------------------------------------------
#READ THE DATA
model_list <- readRDS("fitTargeted_70o.RDS")
test <- readRDS("testSets_70o.RDS")
model_types <- c("glmnet")

#model_types <- c("bridge","blasso","svmLinear","ranger","glmnet")


pdx_cv <- summary_models(model_list = model_list,validation_type = "CV",
               model_types = model_types ,data_source = "PDX_all",test_list = test)


