# import mutation/gene expression/drug response datasets from BeatAML study.
# subset for onkokb features for mutation data, convert gex to gene-set scores
# The script uses input from ./data folder and generates "beatAML.prepared.RDS" file

# https://bioconductor.org/packages/release/bioc/vignettes/PharmacoGx/inst/doc/PharmacoGx.pdf
library(PharmacoGx)
library(data.table)
library(singscore)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = T)

dataDir <- args[1]

# Download BeatAML dataset from orcestra.ca database
# If file already exists, reads from target file
message(date(), " => Importing BeatAML_2018 pSet from orcestra database")
beatAML <- PharmacoGx::downloadPSet(
  name = "BeatAML_2018",
  saveDir = getwd(),
  pSetFileName = "BeatAML.rds"
)

message(date(), " => Reading RNA-Seq data")
# import rna-seq data
gex <- assay(summarizeMolecularProfiles(beatAML,
  cellNames(beatAML),
  mDataType = "rnaseq",
  verbose = FALSE
))
message(date(), " => Reading mutation data")
# import mutation data
mut <- assay(summarizeMolecularProfiles(beatAML,
  cellNames(beatAML),
  mDataType = "mutationall",
  summary.stat = "or",
  verbose = FALSE
))

# check which samples have both gex and mut measurements
to_remove <- unique(c(
  names(which(apply(gex, 2, function(x) sum(is.na(x)) > ncol(gex) / 10))),
  names(which(apply(mut, 2, function(x) sum(is.na(x)) > ncol(mut) / 10)))
))
to_keep <- setdiff(intersect(colnames(gex), colnames(mut)), to_remove)

message(date(), " => Reading drug response data")
# get drug response values
aac <- summarizeSensitivityProfiles(
  beatAML,
  sensitivity.measure = "aac_recomputed",
  summary.stat = "median",
  verbose = FALSE
)

# subset datasets to kept samples
aac <- aac[, to_keep]
mut <- mut[, to_keep]
gex <- gex[, to_keep]


message(date(), " => Converting RNA-seq data to gene-set scores")
# convert gex data to gene-set scores
# 1. remove non protein-coding genes
rd <- rowData(beatAML@molecularProfiles$rnaseq)
gex <- gex[rownames(rd[rd$gene_type == "protein_coding", ]), ]
rownames(gex) <- gsub("\\.[0-9]+", "", rownames(gex))

# get id mapping table from gene ids to gene names
ens2hgnc <- readRDS(file.path(dataDir, "ens2hgnc.RDS"))
# import gene signatures for tumor-microenvironment(TME)-related gene sets from xCell
# and cancer hallmark gene sigatures from MSIGDB
tme <- readRDS(file.path(dataDir, "xCell.signatures.RDS"))
hallmarks <- readRDS(file.path(dataDir, "MSIGDB.hallmarks.RDS"))
geneSets <- c(tme, hallmarks)
geneSets <- sapply(simplify = F, geneSets, function(x) {
  unique(ens2hgnc[match(x, hgnc_symbol)]$ref_gene_id)
})
names(geneSets) <- gsub("\\+", "", names(geneSets))
names(geneSets) <- gsub("-", "_", names(geneSets))
# compute single-sample gene set activity scores
score_gene_set <- function(rankData, genes_up, gene_down = NULL) {
  scores <- singscore::simpleScore(rankData, upSet = genes_up)
  return(scores$TotalScore)
}
rankData <- singscore::rankGenes(gex)
gex_gs <- pbapply::pbsapply(geneSets, function(gs) {
  genes <- intersect(gs, rownames(rankData))
  score_gene_set(rankData, genes)
})
rownames(gex_gs) <- colnames(gex)
gex_gs <- t(gex_gs)

# subset mutation data to only include onkokb genes
onkokb <- readLines(file.path(dataDir, "onkokb_cancerGeneList.txt"))
mut <- mut[intersect(rownames(mut), onkokb), ]
mut <- matrix(as.numeric(mut), ncol = ncol(mut), dimnames = dimnames(mut))

# save datasets (mut+gex_gs)
aac.dt <- data.table(t(aac), keep.rownames = T)
colnames(aac.dt)[1] <- "name"
dat <- list(
  "mut" = mut, "gex_gs" = gex_gs,
  "drugs" = melt.data.table(aac.dt, id.vars = "name")
)

message(date(), " => Saving prepared data")
saveRDS(dat, file = "beatAML.prepared.RDS")

message(date(), " => Finished")
