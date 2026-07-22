#!/usr/bin/env Rscript

# ============================================================
#   scripts/03_qc.R
# Supported methods:
#   MAD, fixed, none
# Optional filters:
#   Ribosomal content, dropout rate, scDblFinder doublets
# Outputs:
#   results/03_qc/_qc.rds
# results/03_qc/QC_summary_.tsv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.msg <- function(verbose = TRUE) {
  function(...) {
    if (isTRUE(verbose)) cat("[INFO]", ..., "\n")
  }
}

.safe_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

.to_logical <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  x %in% c("true", "t", "yes", "y", "1")
}

run_qc <- function(
    x,
    min_feat = 200,
    min_umi = 500,
    mad_n = 5,
    max_mito = 5,
    calc_ribo = FALSE,
    max_ribo = 3,
    calc_drop = FALSE,
    max_drop = 0.95,
    rm_dbl = FALSE,
    dbl_score_thr = 0.5,
    method = c("MAD", "fixed", "none"),
    fixed_thr = list(max_feat = 6000, max_umi = 20000),
    mito_pat = NULL,
    ribo_pat = NULL,
    species = c("human", "mouse"),
    sample_col = "orig.ident",
    outdir = "results/03_qc",
    merge = TRUE,
    parallel = FALSE,
    n.cores = 2,
    verbose = TRUE
) {
  
  log <- .msg(verbose)
  
  method <- match.arg(method)
  species <- match.arg(species)
  
  if (is.null(mito_pat)) {
    mito_pat <- if (species == "human") "^MT-" else "^mt-"
  }
  
  if (is.null(ribo_pat)) {
    ribo_pat <- if (species == "human") "^RP[LS]" else "^Rp[ls]"
  }
  
  .safe_dir(outdir)
  
  if (inherits(x, "Seurat")) {
    obj_list <- list(x)
    single <- TRUE
  } else if (is.list(x) && all(vapply(x, inherits, logical(1), "Seurat"))) {
    obj_list <- x
    single <- FALSE
  } else {
    stop("`x` must be a Seurat object or a list of Seurat objects.", call. = FALSE)
  }
  
  n <- length(obj_list)
  
  log("Starting QC for", n, "sample(s)")
  
  process_one <- function(i) {
    
    t0 <- Sys.time()
    
    obj <- obj_list[[i]]
    
    sample <- if (sample_col %in% colnames(obj@meta.data)) {
      unique(as.character(obj@meta.data[[sample_col]]))[1]
    } else {
      paste0("Sample_", i)
    }
    
    pre_cells <- ncol(obj)
    
    log("Sample:", sample, "| Cells before QC:", pre_cells)
    
    obj$percent_mito <- PercentageFeatureSet(obj, pattern = mito_pat)
    
    if (isTRUE(calc_ribo)) {
      obj$percent_ribo <- PercentageFeatureSet(obj, pattern = ribo_pat)
    }
    
    if (isTRUE(calc_drop)) {
      counts <- GetAssayData(obj, layer = "counts")
      obj$dropout <- Matrix::colSums(counts == 0) / nrow(counts)
      rm(counts)
      invisible(gc())
    }
    
    nfeature_col <- if ("nFeature_RNA" %in% colnames(obj@meta.data)) {
      "nFeature_RNA"
    } else {
      grep("^nFeature_", colnames(obj@meta.data), value = TRUE)[1]
    }
    
    ncount_col <- if ("nCount_RNA" %in% colnames(obj@meta.data)) {
      "nCount_RNA"
    } else {
      grep("^nCount_", colnames(obj@meta.data), value = TRUE)[1]
    }
    
    if (is.na(nfeature_col) || is.na(ncount_col)) {
      stop("Could not find nFeature_* or nCount_* columns.", call. = FALSE)
    }
    
    nf <- obj@meta.data[[nfeature_col]]
    nc <- obj@meta.data[[ncount_col]]
    
    if (method == "MAD") {
      max_feat <- stats::median(nf, na.rm = TRUE) + mad_n * stats::mad(nf, na.rm = TRUE)
      max_umi  <- stats::median(nc, na.rm = TRUE) + mad_n * stats::mad(nc, na.rm = TRUE)
    } else if (method == "fixed") {
      max_feat <- fixed_thr$max_feat %||% Inf
      max_umi  <- fixed_thr$max_umi %||% Inf
    } else {
      max_feat <- Inf
      max_umi  <- Inf
    }
    
    keep <- nf > min_feat &
      nf < max_feat &
      nc > min_umi &
      obj$percent_mito < max_mito
    
    if (isTRUE(calc_ribo)) {
      keep <- keep & obj$percent_ribo < max_ribo
    }
    
    if (isTRUE(calc_drop)) {
      keep <- keep & obj$dropout < max_drop
    }
    
    keep[is.na(keep)] <- FALSE
    
    filt <- obj[, keep]
    
    post_qc_cells <- ncol(filt)
    pct_rm <- round((1 - post_qc_cells / pre_cells) * 100, 2)
    
    log("Retained after QC:", post_qc_cells, "cells | Removed:", pct_rm, "%")
    
    dbl_n <- 0
    
    if (isTRUE(rm_dbl)) {
      
      if (!requireNamespace("scDblFinder", quietly = TRUE)) {
        stop("Package 'scDblFinder' is required when rm_dbl = TRUE.", call. = FALSE)
      }
      
      if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' is required when rm_dbl = TRUE.", call. = FALSE)
      }
      
      log("Running scDblFinder")
      
      sce <- Seurat::as.SingleCellExperiment(filt)
      sce <- scDblFinder::scDblFinder(sce)
      
      filt$scDblFinder_score <- sce$scDblFinder.score
      filt$scDblFinder_class <- sce$scDblFinder.class
      
      dbl_idx <- filt$scDblFinder_score > dbl_score_thr
      dbl_idx[is.na(dbl_idx)] <- FALSE
      
      dbl_n <- sum(dbl_idx)
      
      filt <- filt[, !dbl_idx]
      
      rm(sce)
      invisible(gc())
      
      log("Removed doublets:", dbl_n)
    }
    
    post_final_cells <- ncol(filt)
    
    t1 <- Sys.time()
    
    qc <- data.frame(
      Sample = sample,
      Pre_Cells = pre_cells,
      Post_QC_Cells = post_qc_cells,
      Post_Final_Cells = post_final_cells,
      Removed_QC_Pct = pct_rm,
      Doublets = dbl_n,
      Min_Features = min_feat,
      Min_UMI = min_umi,
      Max_Mito = max_mito,
      Method = method,
      Max_Features = max_feat,
      Max_UMI = max_umi,
      Runtime_sec = round(as.numeric(difftime(t1, t0, units = "secs")), 2),
      stringsAsFactors = FALSE
    )
    
    list(obj = filt, qc = qc)
  }
  
  if (isTRUE(parallel) && n > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' is required for parallel = TRUE.", call. = FALSE)
    }
    
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "obj_list", "sample_col", "min_feat", "min_umi", "mad_n",
        "max_mito", "calc_ribo", "max_ribo", "calc_drop", "max_drop",
        "rm_dbl", "dbl_score_thr", "method", "fixed_thr",
        "mito_pat", "ribo_pat", "%||%"
      ),
      envir = environment()
    )
    
    parallel::clusterEvalQ(cl, {
      library(Seurat)
      library(Matrix)
    })
    
    results <- parallel::parLapply(cl, seq_len(n), process_one)
    
  } else {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    results <- vector("list", n)
    
    for (i in seq_len(n)) {
      results[[i]] <- process_one(i)
      utils::setTxtProgressBar(pb, i)
    }
    
    close(pb)
  }
  
  obj_list <- lapply(results, `[[`, "obj")
  qc_sum <- do.call(rbind, lapply(results, `[[`, "qc"))
  
  param_tag <- paste0(
    "minFeat", min_feat,
    "_minUMI", min_umi,
    "_mito", max_mito,
    "_", method
  )
  
  qc_file <- file.path(outdir, paste0("QC_summary_", param_tag, ".tsv"))
  
  utils::write.table(
    qc_sum,
    qc_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  if (isTRUE(merge) && length(obj_list) > 1) {
    log("Merging QC-filtered objects")
    
    merged <- Reduce(function(a, b) merge(a, b), obj_list)
    
    if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
      merged <- Seurat::JoinLayers(merged)
    }
    
    return(merged)
  }
  
  if (single) {
    return(obj_list[[1]])
  }
  
  obj_list
}

run_qc_step <- function(
    sample,
    rds_in,
    min_feat = 200,
    min_umi = 500,
    max_mito = 5,
    mad_n = 5,
    rm_dbl = FALSE,
    dbl_score = 0.5,
    species = "human",
    outdir = "results/03_qc",
    method = "MAD",
    calc_ribo = FALSE,
    max_ribo = 3,
    calc_drop = FALSE,
    max_drop = 0.95,
    verbose = TRUE
) {
  
  log <- .msg(verbose)
  
  .safe_dir(outdir)
  
  if (!file.exists(rds_in)) {
    stop("Input RDS not found: ", rds_in, call. = FALSE)
  }
  
  log("Reading:", rds_in)
  
  obj <- readRDS(rds_in)
  
  if (!"sample" %in% colnames(obj@meta.data)) {
    obj$sample <- sample
  }
  
  if (!"orig.ident" %in% colnames(obj@meta.data)) {
    obj$orig.ident <- sample
  }
  
  qc_obj <- run_qc(
    x = obj,
    min_feat = min_feat,
    min_umi = min_umi,
    mad_n = mad_n,
    max_mito = max_mito,
    calc_ribo = calc_ribo,
    max_ribo = max_ribo,
    calc_drop = calc_drop,
    max_drop = max_drop,
    rm_dbl = .to_logical(rm_dbl),
    dbl_score_thr = dbl_score,
    method = method,
    species = species,
    sample_col = "sample",
    outdir = outdir,
    merge = FALSE,
    verbose = verbose
  )
  
  rds_out <- file.path(outdir, paste0(sample, "_qc.rds"))
  
  saveRDS(qc_obj, rds_out)
  
  summary_in <- list.files(
    outdir,
    pattern = "^QC_summary_.*\\.tsv$",
    full.names = TRUE
  )
  
  if (length(summary_in) > 0) {
    latest <- summary_in[which.max(file.info(summary_in)$mtime)]
    file.copy(
      latest,
      file.path(outdir, paste0("QC_summary_", sample, ".tsv")),
      overwrite = TRUE
    )
  }
  
  log("Saved:", rds_out)
  
  invisible(qc_obj)
}

if (exists("snakemake")) {
  
  run_qc_step(
    sample    = snakemake@params[["sample"]],
    rds_in    = snakemake@input[["rds"]],
    min_feat  = snakemake@params[["min_feat"]],
    min_umi   = snakemake@params[["min_umi"]],
    max_mito  = snakemake@params[["max_mito"]],
    mad_n     = snakemake@params[["mad_n"]],
    rm_dbl    = snakemake@params[["rm_dbl"]],
    dbl_score = snakemake@params[["dbl_score"]],
    species   = snakemake@params[["species"]],
    outdir    = snakemake@params[["outdir"]] %||% dirname(snakemake@output[["rds"]])
  )
}