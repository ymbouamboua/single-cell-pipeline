#!/usr/bin/env Rscript

# ============================================================
# scripts/04_decontaminate.R
#
# Ambient RNA contamination correction for Seurat objects.
#
# Supported corrected assays:
#   DecontX
#   SoupX
#
# Main output:
#   results/04_clean/<sample>_clean.rds
#
# Compatible with run_pipeline.R:
#
# run_decontam(
#   sample  = s,
#   rds_in  = ".../<sample>_qc.rds",
#   raw_dir = ".../CellRanger/outs",
#   assays  = c("RNA", "DecontX"),
#   outdir  = "results/04_clean"
# )
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.na(x)) return(y)
  x
}

.msg <- function(verbose = TRUE) {
  function(...) {
    if (isTRUE(verbose)) {
      cat("[INFO]", ..., "\n")
    }
  }
}

.warn <- function(...) {
  cat("[WARN]", ..., "\n")
}

.safe_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

.to_logical <- function(x) {
  if (is.logical(x)) return(x)
  
  tolower(as.character(x)) %in%
    c("true", "t", "yes", "y", "1")
}

.is_missing_path <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    is.na(x) ||
    !nzchar(trimws(as.character(x)))
}

.canonical_barcode <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  extracted <- sub(
    "^.*?([ACGTN]+-[0-9]+)$",
    "\\1",
    x,
    perl = TRUE
  )
  
  valid <- grepl(
    "^[ACGTN]+-[0-9]+$",
    extracted
  )
  
  fallback <- sub("^.*[#_]", "", x)
  fallback <- sub("^.*:", "", fallback)
  
  ifelse(valid, extracted, fallback)
}

# -------------------------------------------------------------------------
# Seurat helpers
# -------------------------------------------------------------------------

.get_counts <- function(
    obj,
    assay = "RNA"
) {
  
  if (!assay %in% Assays(obj)) {
    stop(
      "Assay not found in Seurat object: ",
      assay,
      call. = FALSE
    )
  }
  
  Seurat::GetAssayData(
    object = obj,
    assay = assay,
    layer = "counts"
  )
}

.create_assay <- function(counts) {
  
  counts <- methods::as(
    counts,
    "dgCMatrix"
  )
  
  Seurat::CreateAssayObject(
    counts = counts
  )
}

.set_corrected_assay <- function(
    obj,
    counts,
    assay_name,
    normalize = TRUE
) {
  
  counts <- methods::as(
    counts,
    "dgCMatrix"
  )
  
  # Ensure exact feature/cell order.
  counts <- counts[
    rownames(obj),
    colnames(obj),
    drop = FALSE
  ]
  
  obj[[assay_name]] <- .create_assay(counts)
  
  if (isTRUE(normalize)) {
    obj <- Seurat::NormalizeData(
      object = obj,
      assay = assay_name,
      verbose = FALSE
    )
  }
  
  obj
}

# -------------------------------------------------------------------------
# Locate Cell Ranger matrices
# -------------------------------------------------------------------------

.find_matrix <- function(
    cellranger,
    matrix_type = c("raw", "filtered"),
    prefer_h5 = TRUE
) {
  
  matrix_type <- match.arg(matrix_type)
  
  matrix_name <- if (matrix_type == "raw") {
    "raw_feature_bc_matrix"
  } else {
    "filtered_feature_bc_matrix"
  }
  
  h5_name <- paste0(matrix_name, ".h5")
  
  if (
    file.exists(cellranger) &&
    !dir.exists(cellranger) &&
    basename(cellranger) == h5_name
  ) {
    return(
      list(
        type = "h5",
        path = cellranger
      )
    )
  }
  
  if (!dir.exists(cellranger)) {
    return(NULL)
  }
  
  h5 <- list.files(
    cellranger,
    pattern = paste0(
      "^",
      matrix_name,
      "\\.h5$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  
  dirs <- list.dirs(
    cellranger,
    recursive = TRUE,
    full.names = TRUE
  )
  
  dirs <- dirs[
    basename(dirs) == matrix_name
  ]
  
  if (
    isTRUE(prefer_h5) &&
    length(h5)
  ) {
    return(
      list(
        type = "h5",
        path = h5[1]
      )
    )
  }
  
  if (length(dirs)) {
    return(
      list(
        type = "dir",
        path = dirs[1]
      )
    )
  }
  
  if (length(h5)) {
    return(
      list(
        type = "h5",
        path = h5[1]
      )
    )
  }
  
  NULL
}

.read_10x_matrix <- function(input) {
  
  if (is.null(input)) {
    return(NULL)
  }
  
  matrix <- if (input$type == "h5") {
    Seurat::Read10X_h5(
      input$path
    )
  } else {
    Seurat::Read10X(
      input$path
    )
  }
  
  if (is.list(matrix)) {
    
    if ("Gene Expression" %in% names(matrix)) {
      matrix <- matrix[["Gene Expression"]]
    } else {
      matrix <- matrix[[1]]
    }
  }
  
  methods::as(
    matrix,
    "dgCMatrix"
  )
}

.align_raw_to_filtered_genes <- function(
    raw_counts,
    filtered_counts
) {
  
  common_genes <- intersect(
    rownames(filtered_counts),
    rownames(raw_counts)
  )
  
  if (!length(common_genes)) {
    stop(
      "Raw and filtered count matrices have no common genes.",
      call. = FALSE
    )
  }
  
  if (length(common_genes) < nrow(filtered_counts)) {
    .warn(
      nrow(filtered_counts) - length(common_genes),
      "filtered genes are absent from the raw matrix"
    )
  }
  
  list(
    raw = raw_counts[
      common_genes,
      ,
      drop = FALSE
    ],
    filtered = filtered_counts[
      common_genes,
      ,
      drop = FALSE
    ],
    genes = common_genes
  )
}

# -------------------------------------------------------------------------
# DecontX
# -------------------------------------------------------------------------

.run_decontx <- function(
    obj,
    assay = "RNA",
    raw_counts = NULL,
    batch_col = NULL,
    max_iter = 500,
    seed = 1234,
    verbose = TRUE
) {
  
  log <- .msg(verbose)
  
  if (!requireNamespace(
    "SingleCellExperiment",
    quietly = TRUE
  )) {
    stop(
      "Package 'SingleCellExperiment' is required for DecontX.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace(
    "celda",
    quietly = TRUE
  )) {
    stop(
      "Package 'celda' is required for DecontX.",
      call. = FALSE
    )
  }
  
  counts <- .get_counts(
    obj,
    assay = assay
  )
  
  log(
    "Running DecontX on",
    ncol(counts),
    "cells and",
    nrow(counts),
    "genes"
  )
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts
    ),
    colData = obj@meta.data
  )
  
  batch <- NULL
  
  if (
    !is.null(batch_col) &&
    batch_col %in% colnames(obj@meta.data)
  ) {
    batch <- as.factor(
      obj@meta.data[[batch_col]]
    )
  }
  
  background <- NULL
  
  if (!is.null(raw_counts)) {
    
    aligned <- .align_raw_to_filtered_genes(
      raw_counts = raw_counts,
      filtered_counts = counts
    )
    
    raw_aligned <- aligned$raw
    filtered_aligned <- aligned$filtered
    
    # Empty droplets are raw barcodes absent from the filtered object.
    raw_bc <- .canonical_barcode(
      colnames(raw_aligned)
    )
    
    filtered_bc <- .canonical_barcode(
      colnames(filtered_aligned)
    )
    
    ambient_cells <- !raw_bc %in% filtered_bc
    
    if (sum(ambient_cells) >= 50) {
      
      background <- raw_aligned[
        ,
        ambient_cells,
        drop = FALSE
      ]
      
      # Limit extremely large background matrices while retaining
      # representative empty droplets.
      if (ncol(background) > 20000) {
        set.seed(seed)
        
        selected <- sample(
          seq_len(ncol(background)),
          size = 20000
        )
        
        background <- background[
          ,
          selected,
          drop = FALSE
        ]
      }
      
      # DecontX requires matching genes.
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
          counts = filtered_aligned
        ),
        colData = obj@meta.data
      )
      
      log(
        "Using",
        ncol(background),
        "raw droplets as DecontX background"
      )
      
    } else {
      .warn(
        "Fewer than 50 candidate empty droplets were identified.",
        "DecontX will run without an external background."
      )
    }
  }
  
  decontx_args <- list(
    x = sce,
    maxIter = as.integer(max_iter),
    seed = as.integer(seed)
  )
  
  decontx_formals <- names(
    formals(celda::decontX)
  )
  
  if (
    !is.null(batch) &&
    "batch" %in% decontx_formals
  ) {
    decontx_args$batch <- batch
  }
  
  if (
    !is.null(background) &&
    "background" %in% decontx_formals
  ) {
    decontx_args$background <- background
  }
  
  result <- do.call(
    celda::decontX,
    decontx_args
  )
  
  corrected <- SummarizedExperiment::assay(
    result,
    "decontXcounts"
  )
  
  contamination <- SummarizedExperiment::colData(
    result
  )$decontX_contamination
  
  if (is.null(contamination)) {
    contamination <- rep(
      NA_real_,
      ncol(corrected)
    )
  }
  
  list(
    counts = methods::as(
      corrected,
      "dgCMatrix"
    ),
    contamination = contamination,
    sce = result
  )
}

# -------------------------------------------------------------------------
# SoupX
# -------------------------------------------------------------------------

.run_soupx <- function(
    obj,
    raw_counts,
    assay = "RNA",
    cluster_col = "seurat_clusters",
    contamination_fraction = NULL,
    tfidf_min = 1,
    soup_quantile = 0.9,
    force_accept = FALSE,
    verbose = TRUE
) {
  
  log <- .msg(verbose)
  
  if (!requireNamespace(
    "SoupX",
    quietly = TRUE
  )) {
    stop(
      "Package 'SoupX' is required for SoupX correction.",
      call. = FALSE
    )
  }
  
  if (is.null(raw_counts)) {
    stop(
      "Raw Cell Ranger counts are required for SoupX.",
      call. = FALSE
    )
  }
  
  filtered_counts <- .get_counts(
    obj,
    assay = assay
  )
  
  aligned <- .align_raw_to_filtered_genes(
    raw_counts = raw_counts,
    filtered_counts = filtered_counts
  )
  
  raw_counts <- aligned$raw
  filtered_counts <- aligned$filtered
  
  raw_barcode <- .canonical_barcode(
    colnames(raw_counts)
  )
  
  filtered_barcode <- .canonical_barcode(
    colnames(filtered_counts)
  )
  
  raw_match <- match(
    filtered_barcode,
    raw_barcode
  )
  
  if (anyNA(raw_match)) {
    .warn(
      sum(is.na(raw_match)),
      "filtered barcodes were not found in the raw matrix."
    )
  }
  
  # SoupX expects filtered matrix barcodes to correspond to raw barcodes.
  valid_cells <- !is.na(raw_match)
  
  if (sum(valid_cells) < 10) {
    stop(
      "Too few filtered cells matched the raw matrix for SoupX.",
      call. = FALSE
    )
  }
  
  filtered_counts <- filtered_counts[
    ,
    valid_cells,
    drop = FALSE
  ]
  
  filtered_object_cells <- colnames(obj)[valid_cells]
  
  # Use canonical barcode names inside SoupX.
  colnames(filtered_counts) <- filtered_barcode[valid_cells]
  colnames(raw_counts) <- raw_barcode
  
  # Remove duplicated raw barcodes if Cell Ranger paths were combined.
  raw_keep <- !duplicated(colnames(raw_counts))
  
  raw_counts <- raw_counts[
    ,
    raw_keep,
    drop = FALSE
  ]
  
  soup_channel <- SoupX::SoupChannel(
    tod = raw_counts,
    toc = filtered_counts
  )
  
  # Add clustering information.
  if (
    cluster_col %in% colnames(obj@meta.data)
  ) {
    
    clusters <- obj@meta.data[
      filtered_object_cells,
      cluster_col
    ]
    
  } else {
    
    log(
      "Cluster column not found:",
      cluster_col,
      "| Performing temporary clustering for SoupX"
    )
    
    temp <- subset(
      obj,
      cells = filtered_object_cells
    )
    
    temp <- Seurat::NormalizeData(
      temp,
      verbose = FALSE
    )
    
    temp <- Seurat::FindVariableFeatures(
      temp,
      nfeatures = min(
        2000,
        nrow(temp)
      ),
      verbose = FALSE
    )
    
    temp <- Seurat::ScaleData(
      temp,
      verbose = FALSE
    )
    
    temp <- Seurat::RunPCA(
      temp,
      npcs = min(
        30,
        ncol(temp) - 1
      ),
      verbose = FALSE
    )
    
    dims_use <- seq_len(
      min(
        20,
        ncol(Seurat::Embeddings(temp, "pca"))
      )
    )
    
    temp <- Seurat::FindNeighbors(
      temp,
      dims = dims_use,
      verbose = FALSE
    )
    
    temp <- Seurat::FindClusters(
      temp,
      resolution = 0.5,
      verbose = FALSE
    )
    
    clusters <- as.character(
      Seurat::Idents(temp)
    )
    
    rm(temp)
    invisible(gc())
  }
  
  names(clusters) <- filtered_barcode[valid_cells]
  
  soup_channel <- SoupX::setClusters(
    soup_channel,
    clusters
  )
  
  # Supply dimensional coordinates if available.
  if ("umap" %in% Reductions(obj)) {
    
    umap <- Seurat::Embeddings(
      obj,
      "umap"
    )[filtered_object_cells, , drop = FALSE]
    
    rownames(umap) <- filtered_barcode[valid_cells]
    
    soup_channel <- SoupX::setDR(
      soup_channel,
      umap
    )
  }
  
  if (!is.null(contamination_fraction)) {
    
    log(
      "Using fixed SoupX contamination fraction:",
      contamination_fraction
    )
    
    soup_channel <- SoupX::setContaminationFraction(
      soup_channel,
      contaminationFraction =
        as.numeric(contamination_fraction),
      forceAccept = .to_logical(force_accept)
    )
    
  } else {
    
    log("Estimating SoupX contamination fraction")
    
    soup_channel <- SoupX::autoEstCont(
      soup_channel,
      tfidfMin = tfidf_min,
      soupQuantile = soup_quantile,
      forceAccept = .to_logical(force_accept)
    )
  }
  
  corrected <- SoupX::adjustCounts(
    soup_channel,
    roundToInt = TRUE
  )
  
  corrected <- methods::as(
    corrected,
    "dgCMatrix"
  )
  
  # Restore Seurat cell names.
  colnames(corrected) <- filtered_object_cells
  
  # Reconstruct a full object-sized matrix. Cells without raw matches retain
  # their original filtered counts.
  full_corrected <- filtered_counts <- .get_counts(
    obj,
    assay = assay
  )
  
  full_corrected[
    rownames(corrected),
    colnames(corrected)
  ] <- corrected
  
  contamination <- rep(
    NA_real_,
    ncol(obj)
  )
  
  names(contamination) <- colnames(obj)
  
  rho <- soup_channel$metaData$rho
  
  if (!is.null(rho)) {
    
    rho_names <- rownames(
      soup_channel$metaData
    )
    
    rho_match <- match(
      filtered_barcode,
      rho_names
    )
    
    contamination[
      !is.na(rho_match)
    ] <- rho[
      rho_match[!is.na(rho_match)]
    ]
  }
  
  list(
    counts = methods::as(
      full_corrected,
      "dgCMatrix"
    ),
    contamination = contamination,
    soup_channel = soup_channel
  )
}

# -------------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------------

run_decontam <- function(
    sample,
    rds_in,
    raw_dir = NULL,
    assays = c("RNA", "DecontX"),
    input_assay = "RNA",
    preferred_assay = NULL,
    create_clean_alias = TRUE,
    normalize_corrected = TRUE,
    batch_col = NULL,
    cluster_col = "seurat_clusters",
    decontx_max_iter = 500,
    soupx_contamination_fraction = NULL,
    soupx_tfidf_min = 1,
    soupx_soup_quantile = 0.9,
    soupx_force_accept = FALSE,
    seed = 1234,
    outdir = "results/04_clean",
    verbose = TRUE
) {
  
  log <- .msg(verbose)
  
  .safe_dir(outdir)
  
  if (!file.exists(rds_in)) {
    stop(
      "Input RDS not found: ",
      rds_in,
      call. = FALSE
    )
  }
  
  log("===================================")
  log("AMBIENT RNA DECONTAMINATION")
  log("===================================")
  log("Sample:", sample)
  log("Input:", rds_in)
  
  obj <- readRDS(rds_in)
  
  if (!inherits(obj, "Seurat")) {
    stop(
      "Input RDS does not contain a Seurat object.",
      call. = FALSE
    )
  }
  
  if (!input_assay %in% Assays(obj)) {
    stop(
      "Input assay not found: ",
      input_assay,
      call. = FALSE
    )
  }
  
  obj$sample <- sample
  
  assays_requested <- unique(
    toupper(as.character(assays))
  )
  
  supported <- c(
    "RNA",
    "DECONTX",
    "SOUPX"
  )
  
  unsupported <- setdiff(
    assays_requested,
    supported
  )
  
  if (length(unsupported)) {
    stop(
      "Unsupported assay/method(s): ",
      paste(unsupported, collapse = ", "),
      ". Supported values: RNA, DecontX, SoupX.",
      call. = FALSE
    )
  }
  
  initial_counts <- .get_counts(
    obj,
    assay = input_assay
  )
  
  initial_cells <- ncol(initial_counts)
  initial_genes <- nrow(initial_counts)
  initial_umis <- sum(initial_counts)
  
  raw_counts <- NULL
  raw_input <- NULL
  
  needs_raw <- any(
    assays_requested %in% c(
      "DECONTX",
      "SOUPX"
    )
  )
  
  if (
    needs_raw &&
    !.is_missing_path(raw_dir)
  ) {
    
    raw_input <- .find_matrix(
      cellranger = as.character(raw_dir),
      matrix_type = "raw",
      prefer_h5 = TRUE
    )
    
    if (!is.null(raw_input)) {
      
      log(
        "Raw matrix:",
        raw_input$path
      )
      
      raw_counts <- .read_10x_matrix(
        raw_input
      )
      
      log(
        "Raw droplets:",
        ncol(raw_counts),
        "| Raw genes:",
        nrow(raw_counts)
      )
      
    } else {
      .warn(
        "No raw_feature_bc_matrix or raw_feature_bc_matrix.h5 found in:",
        raw_dir
      )
    }
  }
  
  method_results <- list()
  summary_rows <- list()
  
  # -----------------------------------------------------------------------
  # DecontX
  # -----------------------------------------------------------------------
  
  if ("DECONTX" %in% assays_requested) {
    
    decontx_start <- Sys.time()
    
    decontx <- tryCatch(
      .run_decontx(
        obj = obj,
        assay = input_assay,
        raw_counts = raw_counts,
        batch_col = batch_col,
        max_iter = decontx_max_iter,
        seed = seed,
        verbose = verbose
      ),
      error = function(e) {
        .warn(
          "DecontX failed:",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(decontx)) {
      
      obj <- .set_corrected_assay(
        obj = obj,
        counts = decontx$counts,
        assay_name = "DecontX",
        normalize = normalize_corrected
      )
      
      obj$DecontX_contamination <-
        as.numeric(decontx$contamination)
      
      corrected_umis <- sum(
        decontx$counts
      )
      
      summary_rows[["DecontX"]] <- data.frame(
        Sample = sample,
        Method = "DecontX",
        Cells = ncol(decontx$counts),
        Genes = nrow(decontx$counts),
        Original_UMI = initial_umis,
        Corrected_UMI = corrected_umis,
        Removed_UMI = initial_umis - corrected_umis,
        Removed_UMI_Pct = round(
          100 *
            (initial_umis - corrected_umis) /
            initial_umis,
          4
        ),
        Mean_Contamination = mean(
          decontx$contamination,
          na.rm = TRUE
        ),
        Median_Contamination = stats::median(
          decontx$contamination,
          na.rm = TRUE
        ),
        Runtime_sec = round(
          as.numeric(
            difftime(
              Sys.time(),
              decontx_start,
              units = "secs"
            )
          ),
          2
        ),
        stringsAsFactors = FALSE
      )
      
      method_results$DecontX <- decontx
      
      log(
        "Created DecontX assay"
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # SoupX
  # -----------------------------------------------------------------------
  
  if ("SOUPX" %in% assays_requested) {
    
    soupx_start <- Sys.time()
    
    soupx <- tryCatch(
      .run_soupx(
        obj = obj,
        raw_counts = raw_counts,
        assay = input_assay,
        cluster_col = cluster_col,
        contamination_fraction =
          soupx_contamination_fraction,
        tfidf_min = soupx_tfidf_min,
        soup_quantile =
          soupx_soup_quantile,
        force_accept =
          soupx_force_accept,
        verbose = verbose
      ),
      error = function(e) {
        .warn(
          "SoupX failed:",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(soupx)) {
      
      obj <- .set_corrected_assay(
        obj = obj,
        counts = soupx$counts,
        assay_name = "SoupX",
        normalize = normalize_corrected
      )
      
      obj$SoupX_contamination <-
        as.numeric(
          soupx$contamination[
            colnames(obj)
          ]
        )
      
      corrected_umis <- sum(
        soupx$counts
      )
      
      summary_rows[["SoupX"]] <- data.frame(
        Sample = sample,
        Method = "SoupX",
        Cells = ncol(soupx$counts),
        Genes = nrow(soupx$counts),
        Original_UMI = initial_umis,
        Corrected_UMI = corrected_umis,
        Removed_UMI = initial_umis - corrected_umis,
        Removed_UMI_Pct = round(
          100 *
            (initial_umis - corrected_umis) /
            initial_umis,
          4
        ),
        Mean_Contamination = mean(
          soupx$contamination,
          na.rm = TRUE
        ),
        Median_Contamination = stats::median(
          soupx$contamination,
          na.rm = TRUE
        ),
        Runtime_sec = round(
          as.numeric(
            difftime(
              Sys.time(),
              soupx_start,
              units = "secs"
            )
          ),
          2
        ),
        stringsAsFactors = FALSE
      )
      
      method_results$SoupX <- soupx
      
      log(
        "Created SoupX assay"
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # Select preferred corrected assay
  # -----------------------------------------------------------------------
  
  available_corrected <- intersect(
    c("DecontX", "SoupX"),
    Assays(obj)
  )
  
  if (is.null(preferred_assay)) {
    
    preferred_assay <- if (
      "DecontX" %in% available_corrected
    ) {
      "DecontX"
    } else if (
      "SoupX" %in% available_corrected
    ) {
      "SoupX"
    } else {
      input_assay
    }
  }
  
  if (!preferred_assay %in% Assays(obj)) {
    .warn(
      "Preferred assay",
      preferred_assay,
      "is unavailable; using",
      input_assay
    )
    
    preferred_assay <- input_assay
  }
  
  obj$decontamination_method <- preferred_assay
  
  # Optional standardized clean assay.
  if (
    isTRUE(create_clean_alias) &&
    preferred_assay != input_assay
  ) {
    
    preferred_counts <- .get_counts(
      obj,
      assay = preferred_assay
    )
    
    obj <- .set_corrected_assay(
      obj = obj,
      counts = preferred_counts,
      assay_name = "clean",
      normalize = normalize_corrected
    )
    
    DefaultAssay(obj) <- "clean"
    
  } else {
    DefaultAssay(obj) <- preferred_assay
  }
  
  # -----------------------------------------------------------------------
  # Save output
  # -----------------------------------------------------------------------
  
  rds_out <- file.path(
    outdir,
    paste0(sample, "_clean.rds")
  )
  
  saveRDS(
    obj,
    rds_out
  )
  
  if (length(summary_rows)) {
    
    decontam_summary <- do.call(
      rbind,
      summary_rows
    )
    
  } else {
    
    decontam_summary <- data.frame(
      Sample = sample,
      Method = "none",
      Cells = initial_cells,
      Genes = initial_genes,
      Original_UMI = initial_umis,
      Corrected_UMI = initial_umis,
      Removed_UMI = 0,
      Removed_UMI_Pct = 0,
      Mean_Contamination = NA_real_,
      Median_Contamination = NA_real_,
      Runtime_sec = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  decontam_summary$Preferred_Assay <-
    preferred_assay
  
  decontam_summary$Default_Assay <-
    DefaultAssay(obj)
  
  decontam_summary$Raw_Matrix <-
    if (!is.null(raw_input)) {
      raw_input$path
    } else {
      NA_character_
    }
  
  decontam_summary$Output <- rds_out
  
  utils::write.table(
    decontam_summary,
    file.path(
      outdir,
      paste0(
        sample,
        "_decontamination_summary.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # Per-cell contamination estimates.
  contamination_columns <- intersect(
    c(
      "DecontX_contamination",
      "SoupX_contamination"
    ),
    colnames(obj@meta.data)
  )
  
  if (length(contamination_columns)) {
    
    per_cell <- data.frame(
      cell = colnames(obj),
      barcode = .canonical_barcode(
        colnames(obj)
      ),
      sample = obj$sample,
      obj@meta.data[
        ,
        contamination_columns,
        drop = FALSE
      ],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    utils::write.table(
      per_cell,
      file.path(
        outdir,
        paste0(
          sample,
          "_contamination_per_cell.tsv"
        )
      ),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  log(
    "Available assays:",
    paste(
      Assays(obj),
      collapse = ", "
    )
  )
  
  log(
    "Default assay:",
    DefaultAssay(obj)
  )
  
  log(
    "Saved:",
    rds_out
  )
  
  invisible(obj)
}

# -------------------------------------------------------------------------
# Snakemake entry point
# -------------------------------------------------------------------------

if (exists("snakemake")) {
  
  run_decontam(
    sample = snakemake@params[["sample"]],
    rds_in = snakemake@input[["rds"]],
    raw_dir = snakemake@input[["raw_dir"]] %||%
      snakemake@params[["raw_dir"]],
    assays = snakemake@params[["assays"]] %||%
      c("RNA", "DecontX"),
    input_assay =
      snakemake@params[["input_assay"]] %||%
      "RNA",
    preferred_assay =
      snakemake@params[["preferred_assay"]],
    create_clean_alias =
      snakemake@params[["create_clean_alias"]] %||%
      TRUE,
    batch_col =
      snakemake@params[["batch_col"]],
    cluster_col =
      snakemake@params[["cluster_col"]] %||%
      "seurat_clusters",
    outdir = dirname(
      snakemake@output[["rds"]]
    )
  )
}