#!/usr/bin/env Rscript

# ============================================================
#   scripts/01_load.R
# Supported inputs:
#   10X HDF5 or matrix directory
# Filtered or raw feature-barcode matrices
# Output:
#   results/01_loaded/_raw.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

.msg <- function(verbose = TRUE) {
  function(...) if (isTRUE(verbose)) cat("[INFO]", ..., "\n")
}

.safe_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

.find_10x_input <- function(cellranger, use_filtered = TRUE, use_h5 = TRUE) {
  
  matrix_name <- if (isTRUE(use_filtered)) {
    "filtered_feature_bc_matrix"
  } else {
    "raw_feature_bc_matrix"
  }
  
  h5_pattern <- paste0(matrix_name, "\\.h5$")
  
  h5 <- list.files(
    cellranger,
    pattern = h5_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  dirs <- list.dirs(
    cellranger,
    recursive = TRUE,
    full.names = TRUE
  )
  
  dirs <- dirs[basename(dirs) == matrix_name]
  
  if (isTRUE(use_h5) && length(h5) > 0) {
    return(list(type = "h5", path = h5[1]))
  }
  
  if (length(dirs) > 0) {
    return(list(type = "dir", path = dirs[1]))
  }
  
  if (length(h5) > 0) {
    return(list(type = "h5", path = h5[1]))
  }
  
  stop("No valid 10X input found in: ", cellranger, call. = FALSE)
}

read_10x_counts <- function(cellranger, use_filtered = TRUE, use_h5 = TRUE) {
  
  input <- .find_10x_input(
    cellranger = cellranger,
    use_filtered = use_filtered,
    use_h5 = use_h5
  )
  
  if (input$type == "h5") {
    counts <- Seurat::Read10X_h5(input$path)
  } else {
    counts <- Seurat::Read10X(input$path)
  }
  
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  
  counts
}

run_load <- function(
    sample,
    cellranger,
    species = "human",
    use_filtered = TRUE,
    use_h5 = TRUE,
    outdir = "results/01_loaded",
    min_cells = 3,
    min_features = 0,
    verbose = TRUE
) {
  
  log <- .msg(verbose)
  .safe_dir(outdir)
  
  if (!dir.exists(cellranger)) {
    stop("Cell Ranger directory not found: ", cellranger, call. = FALSE)
  }
  
  log("Loading sample:", sample)
  log("Input:", cellranger)
  
  counts <- read_10x_counts(
    cellranger = cellranger,
    use_filtered = use_filtered,
    use_h5 = use_h5
  )
  
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = min_cells,
    min.features = min_features
  )
  
  obj$sample <- sample
  obj$orig.ident <- sample
  obj$species <- species
  
  out <- file.path(outdir, paste0(sample, "_raw.rds"))
  
  saveRDS(obj, out)
  
  stats <- data.frame(
    Sample = sample,
    Cells = ncol(obj),
    Genes = nrow(obj),
    Species = species,
    Input = cellranger,
    Output = out,
    stringsAsFactors = FALSE
  )
  
  write.table(
    stats,
    file.path(outdir, paste0(sample, "_load_summary.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  log("Saved:", out)
  invisible(obj)
}

if (exists("snakemake")) {
  run_load(
    sample       = snakemake@params[["sample"]],
    cellranger   = snakemake@input[["cellranger"]],
    species      = snakemake@params[["species"]],
    use_filtered = snakemake@params[["use_filt"]],
    outdir       = dirname(snakemake@output[["rds"]])
  )
}