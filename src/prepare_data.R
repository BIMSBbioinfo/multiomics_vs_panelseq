#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)


# {SETUP}
## libraries & functions
library(tictoc)
library(data.table)
utility_script <- '/data/local/buyar/arcas/uyar_et_al_multiomics_deeplearning/src/common_functions.R' # utility script
source(utility_script)
## time functions
my.msg.tic <- function(tic, msg)
{
  if (is.null(msg) || is.na(msg) || length(msg) == 0)
  {
    outmsg <- paste0("Finished\nElapsed time: ", lubridate::seconds_to_period(round(toc - tic, 0)))
  }
  else
  {
    outmsg <- paste0("Starting ", msg, "...")
  }
}
my.msg.toc <- function(tic, toc, msg, info)
{
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

## function to get gene set scores from gex data 
score_gene_set <- function(rankData, genes_up, gene_down = NULL) {
  scores <- singscore::simpleScore(rankData, upSet = genes_up)
  return(scores$TotalScore)
}



## Paths
dset <- as.character(args[1])
p.data <- paste0("/local/abarano/Projects/DrugResponse/data/Raw/", dset)
p.aux <- "/local/abarano/Projects/DrugResponse/data/auxiliary"
p.out <- "/local/abarano/Projects/DrugResponse/data/"

# {MAIN}
# time stamp
tic(msg = "Gene expression", quiet = FALSE, func.tic = my.msg.tic)
## Load gex data
gex <- data.table::fread(file.path(p.data, 'gex.tsv'))
gex <- as.matrix(data.frame(gex[,-1], row.names = gex$gene_id, check.names = F))
## Load gene sets
### TME-related genesets
tme <- read_xCell_geneSets(file.path(p.aux, "gene_lists", "xCell_cell_signature_genes.tsv"), 
                           group_by_cell_type = T)
### Cancer Hallmarks gene sets
hallmarks <- readRDS(file.path(p.aux, "gene_lists", "hallmark_genesets.all.v7.0.symbols.RDS"))
### Bind & map symbols to ids
geneSets <- c(tme, hallmarks)
ens2hgnc <- readRDS(file.path(p.aux, "ens2hgnc.RDS"))
geneSets <- sapply(simplify = F, geneSets, function(x) {
  unique(ens2hgnc[match(x, hgnc_symbol)]$ref_gene_id)
})
names(geneSets) <- gsub("\\+", "", names(geneSets))
names(geneSets) <- gsub("-", "_", names(geneSets))
## Score gex data
rankData <- singscore::rankGenes(gex)
gex_gs <- pbapply::pbsapply(geneSets, function(gs) {
  genes <- intersect(gs, rownames(rankData))
  score_gene_set(rankData, genes)
})
rownames(gex_gs) <- colnames(gex)
gex_gs <- t(gex_gs)
# time stamp
toc(quiet = FALSE, func.toc = my.msg.toc)


## Load CNV & MUT data, subset for oncokb genes, and summarize per gene
### load & prepate oncokb panel genes
oncokb <- readLines(file.path(p.aux, "gene_lists", "oncokb_cancerGeneList.txt"))
oncokb <- unique(ens2hgnc[match(oncokb, hgnc_symbol)]$ref_gene_id)
oncokb <- oncokb[!is.na(oncokb)]

# time stamp
tic(msg = "Mutation data", quiet = FALSE, func.tic = my.msg.tic)
### load MUT
mut <- data.table::fread(file.path(p.data, 'mut.tsv'))
#### count mutations per gene per sample
if (dset == "PDX") {
  # This particular table kills data.table when i try to use melt -- something really weird
  #dt <- mut[!is.na(Gene)][!is.na(HGVSp)][, length(Gene), by = c('Gene', 'Sample')]
  #dtc <- dcast.data.table(dt, Gene ~ Sample, value.var = 'V1', fill = 0)
  #setnames(dtc, old = "Gene", new = "gene_id")
  library(tidyverse)
  dt <- mut %>% 
    dplyr::filter(!is.na(Gene), !is.na(HGVSp)) %>% 
    dplyr::group_by(Sample, Gene) %>% 
    dplyr::summarize(V1 = n(), .groups = "drop")
  dtc <- dt %>% 
    pivot_wider(., names_from = Sample, values_from = V1, values_fill = 0) %>% 
    dplyr::rename(gene_id = Gene)
} else {
  dt <- mut[!is.na(gene_id)][Variant_Classification != 'Silent'][, length(Chromosome), by = c('gene_id', 'Tumor_Sample_Barcode')]
  dtc <- dcast.data.table(dt, gene_id ~ Tumor_Sample_Barcode, value.var = 'V1', fill = 0)
}
dfc <- data.frame(dtc[,-1], row.names = dtc$gene_id, check.names = F)
mut <- as.matrix(dfc)
mut.panel <- as.matrix(dfc[intersect(rownames(dfc), oncokb), ])

# time stamp
toc(quiet = FALSE, func.toc = my.msg.toc)

# time stamp
tic(msg = "CNV data", quiet = FALSE, func.tic = my.msg.tic)
### load CNV
cnv <- data.table::fread(file.path(p.data, 'cnv.tsv'))
cnv <- as.matrix(data.frame(cnv[,-1], row.names = cnv$gene_id, check.names = F))
cnv.panel <- cnv[intersect(rownames(cnv), oncokb),]
# time stamp
toc(quiet = FALSE, func.toc = my.msg.toc)


# time stamp
tic(msg = "Writing results", quiet = FALSE, func.tic = my.msg.tic)
## Intersect samples between across three data entities & write the output
selected <- intersect(colnames(cnv), intersect(colnames(gex_gs), colnames(mut)))
dat <- list('mut' = mut[, selected], 
            'mut.panel' = mut.panel[, selected],
            'cnv' = cnv[, selected],
            'cnv.panel' = cnv.panel[, selected],
            'gex' = gex_gs[, selected])
saveRDS(dat, file = file.path(p.out, paste0("data_", dset, ".RDS")))
# time stamp
toc(quiet = FALSE, func.toc = my.msg.toc)


