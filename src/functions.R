loadCNV <- function(path_data, dataset, group) {
  # prepare matrix
  tmp <- data.table::fread(file.path(path_data, dataset, group, "cnv.tsv.gz"), data.table = F)
  m <- as.matrix(tmp[, -1])
  rownames(m) <- tmp[["V1"]]

  # get some stats
  v.mns <- matrixStats::colMeans2(m)
  v.vrs <- matrixStats::colVars(m)

  # assemble output
  l.out <- list()
  l.out[["data_raw"]] <- m
  l.out[["means"]] <- v.mns
  l.out[["variances"]] <- v.vrs
  l.out[["diagnostics"]] <- ggplot(tibble(mean = v.mns, variance = v.vrs), aes(mean, variance)) +
    geom_hex() +
    theme_bw(base_size = 13) +
    theme(aspect.ratio = 1)
  attr(l.out, "dataset") <- dataset
  attr(l.out, "group") <- group

  # clean after yourself
  rm(tmp, m, v.mns, v.vrs)
  gc()

  return(l.out)
}


filtCNV <- function(l.dat, perc_cutoff = 0.9) {
  # filter the data
  l.dat[["data_varFiltered"]] <- l.dat[["data_raw"]][, which(l.dat[["variances"]] >= quantile(l.dat[["variances"]], perc_cutoff)) + 1]
  attr(l.dat, "variance_cutoff") <- quantile(l.dat[["variances"]], perc_cutoff)
  return(l.dat)
}


loadMut <- function(path_data, dataset, group) {
  # little data processing
  tmp <- data.table::fread(file.path(path_data, dataset, group, "mut.tsv.gz"), data.table = F)
  if (!str_detect(str_to_lower(dataset), "pdx")) {
    ayy <- tmp %>%
      rename(cell_lines = sample_id) %>%
      select(cell_lines, gene_id, protein_change) %>%
      filter(!protein_change == "p.?") %>%
      distinct() %>%
      reshape::melt(
        id.vars = c("cell_lines", "gene_id"),
        na.rm = TRUE, nvalue.name = "protein_change"
      ) %>%
      select(!variable) %>%
      filter(value != "") %>%
      mutate(value = gsub("p.", "", value)) %>%
      pivot_wider(.,
        names_from = c("gene_id", "value"),
        values_from = value, values_fn = length, values_fill = 0
      )
  } else {
    ayy <- tmp %>%
      rename(gene_id = Gene, cell_lines = sample_id) %>%
      select(cell_lines, gene_id, protein_change) %>%
      filter(
        !protein_change == "p.?",
        !str_starts(protein_change, "-")
      ) %>%
      distinct() %>%
      reshape::melt(
        id.vars = c("cell_lines", "gene_id"),
        na.rm = TRUE, nvalue.name = "protein_change"
      ) %>%
      select(!variable) %>%
      filter(value != "") %>%
      mutate(value = gsub("p.", "", value)) %>%
      pivot_wider(.,
        names_from = c("gene_id", "value"),
        values_from = value, values_fn = length, values_fill = 0
      )
  }
  # assmeble a matrix
  m <- as.matrix(ayy[, -1])
  rownames(m) <- ayy[["cell_lines"]]

  # get some stats
  v.mns <- matrixStats::colMeans2(m)
  v.vrs <- matrixStats::colVars(m)

  # assemble output
  l.out <- list()
  l.out[["data_raw"]] <- m
  l.out[["means"]] <- v.mns
  l.out[["variances"]] <- v.vrs
  l.out[["diagnostics"]] <- ggplot(tibble(mean = v.mns, variance = v.vrs), aes(mean, variance)) +
    geom_hex() +
    theme_bw(base_size = 13) +
    theme(aspect.ratio = 1)
  attr(l.out, "dataset") <- dataset
  attr(l.out, "group") <- group

  # clean after yourself
  rm(tmp, ayy, m, v.mns, v.vrs)
  gc()

  return(l.out)
}


sumMutPerGene <- function(l.dat) {
  if (dim(l.dat$data_raw)[2] > 3000) {
    ## split data into chunks and summarize one by one
    ### select number of chunks
    itrtr <- seq(2, 100, by = 2)
    v.sel <- abs((dim(l.dat$data_raw)[2] / itrtr) - 1500)
    divr <- itrtr[which(v.sel == min(v.sel))]
    
    ## iterate
    tmp <- tibble()
    for (i in 1:divr) {
      l.dat$data_raw[, (1500 * (i - 1) + 1):ifelse(i == divr, dim(l.dat$data_raw)[2], 1500 * i)] %>% 
        as_tibble(., rownames = "sample_id") %>% 
        pivot_longer(., contains("ENS"), names_to = "feature", values_to = "value") %>% 
        filter(value > 0) %>% 
        mutate(feature = str_split_fixed(feature, "_", 2)[, 1]) %>% 
        group_by(sample_id, feature) %>% 
        summarize(mut_count = sum(value), .groups = "drop") %>% 
        #pivot_wider(., names_from = "feature", values_from = "mut_count") %>% 
        rbind(tmp, .) -> tmp
      #
      tmp <- tmp %>% 
        group_by(sample_id, feature) %>% 
        summarize(mut_count = sum(mut_count), .groups = "drop")
    }
    # spread
    tmp <- tmp %>% pivot_wider(., names_from = "feature", values_from = "mut_count", values_fill = 0)
    # turn back into matrix
    rnms <- tmp$sample_id
    tmp <- as.matrix(tmp[, -1])
    rownames(tmp) <- rnms
    # update list
    l.dat$data_raw <- tmp
  } else {
    l.dat$data_raw %>% 
      as_tibble(., rownames = "sample_id") %>% 
      pivot_longer(., contains("ENS"), names_to = "feature", values_to = "value") %>% 
      filter(value > 0) %>% 
      mutate(feature = str_split_fixed(feature, "_", 2)[, 1]) %>% 
      group_by(sample_id, feature) %>% 
      summarize(mut_count = sum(value), .groups = "drop") %>% 
      pivot_wider(., names_from = "feature", values_from = "mut_count", values_fill = 0) -> tmp
    # turn back into matrix
    rnms <- tmp$sample_id
    tmp <- as.matrix(tmp[, -1])
    rownames(tmp) <- rnms
    # update list
    l.dat$data_raw <- tmp
  }
  return(l.dat)
}



filtMut <- function(l.dat, var_coef = -0.5, var_intercept = 0.05) {
  # filter the data
  l.dat[["data_varFiltered"]] <- l.dat[["data_raw"]][, which(l.dat[["variances"]] >= var_coef * l.dat[["means"]] + var_intercept) + 1]
  attr(l.dat, "variance_cutoff") <- paste0("Variance >= ", var_coef, "* Mean + ", var_intercept)
  return(l.dat)
}



loadGEX <- function(path_data, dataset, group,
                    m.path_depmap = "",
                    m.path_pdx = "data/PDX_all/meta_info.tsv") {
  # prepare matrix
  tmp <- data.table::fread(file.path(path_data, dataset, group, "gex.tsv.gz"),
    data.table = F
  )
  m <- as.matrix(tmp[, -1])
  rownames(m) <- tmp[["V1"]]

  if (str_detect(str_to_lower(dataset), "pdx")) {
    m <- apply(m, 2, function(x) log2(x + 1))
  }
  
  # get some stats
  v.mns <- matrixStats::colMeans2(m)
  v.vrs <- matrixStats::colVars(m)
  
  ## filter by baseline expression -> > 5 log2(expr + 1)
  m <- m[, which(v.mns > 5)]
  
  # load meta data for qsmooth
  if (!str_detect(str_to_lower(dataset), "pdx")) {
    t.m <- data.table::fread("/local/abarano/Projects/DrugResponse/data/auxiliary/CCLE/depmap_sample_info.csv",
                             data.table = F
    ) %>%
      dplyr::select(sample_id = DepMap_ID, tumor_type = primary_disease)
  } else {
    t.m <- data.table::fread("/local/abarano/Projects/DrugResponse/data/auxiliary/PDX/meta_info.tsv",
                             data.table = F
    ) %>%
      filter(column_name == "tumor_type") %>%
      dplyr::select(sample_id, tumor_type = value) %>%
      distinct(sample_id, .keep_all = T)
  }
  
  # transpose the matrix
  m.t <- t(m)
  
  # filter meta
  samp.int <- intersect(colnames(m.t), t.m$sample_id)
  # match the samples
  t.m <- t.m[t.m$sample_id %in% samp.int, ]
  m.t <- m.t[, samp.int]
  
  # normalize with qsmooth for future use
  t.m_norm <- qsmooth(object = m.t, group_factor = t.m$tumor_type)
  m <- t(as.matrix(t.m_norm@qsmoothData))
  
  # assemble output
  l.out <- list()
  l.out[["data_raw"]] <- m
  l.out[["means"]] <- v.mns
  l.out[["variances"]] <- v.vrs
  l.out[["diagnostics"]] <- ggplot(
    tibble(mean = v.mns, variance = v.vrs),
    aes(mean, variance)
  ) +
    geom_hex() +
    geom_vline(xintercept = 5, color = "red") +
    theme_bw(base_size = 13) +
    theme(aspect.ratio = 1)
  attr(l.out, "dataset") <- dataset
  attr(l.out, "group") <- group

  # clean after yourself
  rm(tmp, m, v.mns, v.vrs)
  gc()

  return(l.out)
}




filtGEX <- function(l.dat,
                    var_coef = -1.2,
                    var_intercept = 7.5,
                    perc_cutoff = 0.8) {
  v.vrs <- l.dat$variances
  v.mns <- l.dat$means
  # filter the data by variance
  l.dat[["data_varFiltered"]] <- m[, which(v.vrs >= var_coef * v.mns + var_intercept & v.vrs >= quantile(v.vrs, 0.9)) + 1]
  attr(l.dat, "variance_cutoff") <- paste0(
    "Variance >= ", var_coef, "* Mean + ",
    var_intercept, "; ", perc_cutoff
  )
  l.dat[["diagnostic_cutoffs"]] <- ggplot(
    tibble(mean = v.mns, variance = v.vrs),
    aes(mean, variance)
  ) +
    geom_hex() +
    geom_abline(intercept = var_intercept, slope = var_coef) +
    geom_hline(yintercept = quantile(v.vrs, perc_cutoff)) +
    theme_bw(base_size = 13) +
    theme(aspect.ratio = 1)

  # clean after yourself
  rm(m.t, t.m, samp.int, t.m_norm, v.mns, v.vrs)
  gc()
  return(l.dat)
}
