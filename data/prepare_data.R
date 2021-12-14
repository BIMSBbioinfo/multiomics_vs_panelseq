# Remove everything else
rm(list = ls())
gc()

# This script is written to prepare drug response data

library(tidyr)
library(dplyr)

# paths
p.d <- "data"
l.d <- c("CCLE", "PDX_all")

# 1. read drug response files
# 2. Find each drug
# 3. Remove missing values and sample_id
# 4. Extract only unique sample_ids

for (source in l.d) {
  responses <- furrr::future_map(
    list.files(paste0(p.d, "/", source, "/raw_drug_response"),
      pattern = "*.tsv", full.names = TRUE
    ),
    ~ data.table::fread(.x, data.table = F)
  ) %>%
    `names<-`(source)
  drugs <- unique(responses[[source]]$column_name)
  readr::write_tsv(
    as.data.frame(t(drugs)),
    paste0("misc", "/drugs_info", ".tsv")
  )
  for (drug in drugs) {
    res <- responses[[source]] %>%
      filter(column_name == drug) %>%
      filter(!is.na(value) & !is.na(sample_id))
    res <- res[res$sample_id %in% unique(res$sample_id), ]
    readr::write_tsv(res, paste0(p.d, "/", source, "/drug_response/", drug, ".tsv"))
  }
}
