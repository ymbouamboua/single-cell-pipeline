#!/usr/bin/env Rscript

# ============================================================
# scripts/05_integrate.R
#
# Supported methods:
#   MERGE / NONE, HARMONY, RPCA, CCA
#
# Output:
#   results/05_integrated/integrated.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

msg  <- function(...) cat("[INFO]", ..., "\n")
warn <- function(...) cat("[WARN]", ..., "\n")

.safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.clean_sample_name <- function(path) {
  x <- tools::file_path_sans_ext(basename(path))
  sub("_(clean|qc|demuxed|raw)$", "", x)
}

.prepare_object <- function(object, sample, batch_col) {
  if (!inherits(object, "Seurat")) {
    stop("RDS is not a Seurat object: ", sample, call. = FALSE)
  }
  
  if (!ncol(object)) {
    stop("Seurat object contains no cells: ", sample, call. = FALSE)
  }
  
  if (!"sample" %in% colnames(object@meta.data)) {
    object$sample <- sample
  }
  
  if (!"orig.ident" %in% colnames(object@meta.data)) {
    object$orig.ident <- sample
  }
  
  if (!batch_col %in% colnames(object@meta.data)) {
    stop(
      "Missing batch column ", dQuote(batch_col),
      " in ", dQuote(sample), ". Available columns: ",
      paste(colnames(object@meta.data), collapse = ", "),
      call. = FALSE
    )
  }
  
  batch <- trimws(as.character(object@meta.data[[batch_col]]))
  
  invalid <- is.na(batch) |
    !nzchar(batch) |
    tolower(batch) %in% c(
      "na", "nan", "null", "none",
      "unknown", "unassigned"
    )
  
  if (any(invalid)) {
    stop(
      "Batch column ", dQuote(batch_col),
      " contains ", sum(invalid),
      " invalid value(s) in ", dQuote(sample), ".",
      call. = FALSE
    )
  }
  
  object@meta.data[[batch_col]] <- factor(batch)
  object
}

.merge_objects <- function(objects, sample_names) {
  if (length(objects) == 1L) {
    return(objects[[1L]])
  }
  
  merge(
    x = objects[[1L]],
    y = objects[-1L],
    add.cell.ids = sample_names,
    project = "scRNA_integrated"
  )
}

.join_layers <- function(object) {
  if (
    exists(
      "JoinLayers",
      envir = asNamespace("SeuratObject"),
      mode = "function"
    )
  ) {
    object <- SeuratObject::JoinLayers(object)
  }
  
  object
}

.validate_dims <- function(dims, npcs) {
  dims <- as.integer(dims[[1L]])
  npcs <- as.integer(npcs[[1L]])
  
  if (is.na(dims) || dims < 1L) {
    stop("'dims' must be a positive integer.", call. = FALSE)
  }
  
  if (is.na(npcs) || npcs < 2L) {
    stop("'npcs' must be greater than 1.", call. = FALSE)
  }
  
  if (dims > npcs) {
    warn("'dims' exceeds 'npcs'; using", npcs, "dimensions.")
    dims <- npcs
  }
  
  seq_len(dims)
}

.preprocess_object <- function(
    object,
    nfeatures,
    npcs,
    seed,
    features = NULL
) {
  object <- NormalizeData(object, verbose = FALSE)
  
  object <- FindVariableFeatures(
    object,
    nfeatures = nfeatures,
    verbose = FALSE
  )
  
  object <- ScaleData(
    object,
    features = features,
    verbose = FALSE
  )
  
  RunPCA(
    object,
    features = features,
    npcs = npcs,
    seed.use = seed,
    verbose = FALSE
  )
}

.cluster_object <- function(
    object,
    reduction,
    dims,
    resolution,
    seed
) {
  object <- FindNeighbors(
    object,
    reduction = reduction,
    dims = dims,
    verbose = FALSE
  )
  
  object <- FindClusters(
    object,
    resolution = resolution,
    random.seed = seed,
    verbose = FALSE
  )
  
  RunUMAP(
    object,
    reduction = reduction,
    dims = dims,
    seed.use = seed,
    verbose = FALSE
  )
}

.write_tsv <- function(x, file, row.names = FALSE, col.names = TRUE) {
  utils::write.table(
    x,
    file = file,
    sep = "\t",
    row.names = row.names,
    col.names = col.names,
    quote = FALSE
  )
}

.save_outputs <- function(
    object,
    outdir,
    method,
    batch_col,
    dims_use,
    resolution,
    nfeatures,
    npcs,
    seed,
    n_objects
) {
  integrated_file <- file.path(outdir, "integrated.rds")
  saveRDS(object, integrated_file)
  
  plot_group <- if (batch_col %in% colnames(object@meta.data)) {
    batch_col
  } else {
    "sample"
  }
  
  p1 <- DimPlot(
    object,
    reduction = "umap",
    group.by = plot_group,
    raster = TRUE
  ) +
    ggtitle(paste("UMAP by", plot_group))
  
  p2 <- DimPlot(
    object,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("UMAP by cluster")
  
  grDevices::pdf(
    file.path(outdir, "umap_overview.pdf"),
    width = 12,
    height = 6
  )
  
  print(p1 + p2)
  grDevices::dev.off()
  
  cluster_summary <- as.data.frame(
    table(cluster = Idents(object)),
    stringsAsFactors = FALSE
  )
  
  colnames(cluster_summary) <- c("cluster", "cells")
  
  .write_tsv(
    cluster_summary,
    file.path(outdir, "cluster_summary.tsv")
  )
  
  sample_col <- if ("sample" %in% colnames(object@meta.data)) {
    "sample"
  } else {
    batch_col
  }
  
  cluster_by_sample <- as.data.frame.matrix(
    table(
      object@meta.data[[sample_col]],
      Idents(object)
    )
  )
  
  .write_tsv(
    cluster_by_sample,
    file.path(outdir, "cluster_by_sample.tsv"),
    row.names = TRUE,
    col.names = NA
  )
  
  parameters <- data.frame(
    parameter = c(
      "method", "batch_col", "dims",
      "resolution", "nfeatures", "npcs",
      "seed", "objects", "cells"
    ),
    value = c(
      method, batch_col, length(dims_use),
      resolution, nfeatures, npcs,
      seed, n_objects, ncol(object)
    ),
    stringsAsFactors = FALSE
  )
  
  .write_tsv(
    parameters,
    file.path(outdir, "integration_parameters.tsv")
  )
  
  msg("Integrated cells:", ncol(object))
  msg("Integrated features:", nrow(object))
  msg("Saved:", integrated_file)
  
  invisible(object)
}

run_integration <- function(
    rds_paths,
    method = "Harmony",
    batch_col = "status",
    dims = 30,
    resolution = 0.5,
    nfeatures = 3000,
    npcs = 50,
    seed = 1234,
    outdir = "results/05_integrated",
    verbose = TRUE
) {
  outdir <- .safe_dir(outdir)
  
  rds_paths <- trimws(as.character(rds_paths))
  rds_paths <- rds_paths[!is.na(rds_paths) & nzchar(rds_paths)]
  
  if (!length(rds_paths)) {
    stop("No RDS files provided.", call. = FALSE)
  }
  
  missing <- rds_paths[!file.exists(rds_paths)]
  
  if (length(missing)) {
    stop(
      "Missing RDS file(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  method <- toupper(trimws(as.character(method[[1L]])))
  batch_col <- trimws(as.character(batch_col[[1L]]))
  
  if (!nzchar(batch_col)) {
    batch_col <- "status"
  }
  
  resolution <- as.numeric(resolution[[1L]])
  nfeatures <- as.integer(nfeatures[[1L]])
  npcs <- as.integer(npcs[[1L]])
  seed <- as.integer(seed[[1L]])
  
  if (is.na(seed)) {
    seed <- 1234L
  }
  
  if (is.na(resolution) || resolution <= 0) {
    stop("'resolution' must be positive.", call. = FALSE)
  }
  
  if (is.na(nfeatures) || nfeatures < 1L) {
    stop("'nfeatures' must be positive.", call. = FALSE)
  }
  
  dims_use <- .validate_dims(dims, npcs)
  set.seed(seed)
  
  sample_names <- vapply(
    rds_paths,
    .clean_sample_name,
    character(1L)
  )
  
  duplicated_samples <- unique(
    sample_names[duplicated(sample_names)]
  )
  
  if (length(duplicated_samples)) {
    stop(
      "Duplicated sample names: ",
      paste(duplicated_samples, collapse = ", "),
      call. = FALSE
    )
  }
  
  objects <- lapply(rds_paths, readRDS)
  names(objects) <- sample_names
  
  objects <- Map(
    function(object, sample) {
      .prepare_object(object, sample, batch_col)
    },
    objects,
    sample_names
  )
  
  batch_levels <- sort(unique(unlist(
    lapply(
      objects,
      function(x) as.character(x@meta.data[[batch_col]])
    ),
    use.names = FALSE
  )))
  
  if (isTRUE(verbose)) {
    msg("===================================")
    msg("SCRNA INTEGRATION")
    msg("===================================")
    msg("Method:", method)
    msg("Batch column:", batch_col)
    msg("Input objects:", length(objects))
    msg("Dimensions:", paste(range(dims_use), collapse = "-"))
    msg("Resolution:", resolution)
    msg("Batch levels:", paste(batch_levels, collapse = ", "))
  }
  
  stats <- do.call(
    rbind,
    Map(
      function(object, sample) {
        counts <- table(
          object@meta.data[[batch_col]],
          useNA = "ifany"
        )
        
        data.frame(
          sample = sample,
          batch_col = batch_col,
          batch = names(counts),
          cells = as.integer(counts),
          genes = nrow(object),
          stringsAsFactors = FALSE
        )
      },
      objects,
      sample_names
    )
  )
  
  rownames(stats) <- NULL
  
  .write_tsv(
    stats,
    file.path(outdir, "integration_stats.tsv")
  )
  
  if (length(objects) == 1L) {
    warn("Only one object supplied; using MERGE.")
    method <- "MERGE"
  }
  
  if (method %in% c("NONE", "MERGE")) {
    msg("Merging without batch correction")
    
    object <- .merge_objects(objects, sample_names)
    
    object <- .preprocess_object(
      object,
      nfeatures,
      npcs,
      seed
    )
    
    object <- .cluster_object(
      object,
      reduction = "pca",
      dims = dims_use,
      resolution = resolution,
      seed = seed
    )
    
  } else if (method %in% c("RPCA", "CCA")) {
    msg("Running Seurat integration:", method)
    
    objects <- lapply(
      objects,
      function(x) {
        x <- NormalizeData(x, verbose = FALSE)
        
        FindVariableFeatures(
          x,
          nfeatures = nfeatures,
          verbose = FALSE
        )
      }
    )
    
    features <- SelectIntegrationFeatures(
      object.list = objects,
      nfeatures = nfeatures
    )
    
    objects <- lapply(
      objects,
      .preprocess_object,
      nfeatures = nfeatures,
      npcs = npcs,
      seed = seed,
      features = features
    )
    
    anchors <- FindIntegrationAnchors(
      object.list = objects,
      anchor.features = features,
      reduction = tolower(method),
      dims = dims_use
    )
    
    object <- IntegrateData(
      anchorset = anchors,
      dims = dims_use
    )
    
    DefaultAssay(object) <- "integrated"
    
    object <- ScaleData(
      object,
      verbose = FALSE
    )
    
    object <- RunPCA(
      object,
      npcs = npcs,
      seed.use = seed,
      verbose = FALSE
    )
    
    object <- .cluster_object(
      object,
      reduction = "pca",
      dims = dims_use,
      resolution = resolution,
      seed = seed
    )
    
  } else if (method == "HARMONY") {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop(
        "Package 'harmony' is required. Install it with:\n",
        "install.packages('harmony')",
        call. = FALSE
      )
    }
    
    if (length(batch_levels) < 2L) {
      stop(
        "Harmony requires at least two levels in ",
        dQuote(batch_col), ".",
        call. = FALSE
      )
    }
    
    msg("Running Harmony using:", batch_col)
    
    object <- .merge_objects(objects, sample_names)
    object <- .join_layers(object)
    
    object <- .preprocess_object(
      object,
      nfeatures,
      npcs,
      seed
    )
    
    object <- harmony::RunHarmony(
      object = object,
      group.by.vars = batch_col,
      reduction.use = "pca",
      dims.use = dims_use,
      verbose = FALSE
    )
    
    object <- .cluster_object(
      object,
      reduction = "harmony",
      dims = dims_use,
      resolution = resolution,
      seed = seed
    )
    
  } else {
    stop(
      "Unsupported method: ", method,
      ". Use MERGE, NONE, HARMONY, RPCA, or CCA.",
      call. = FALSE
    )
  }
  
  if ("RNA" %in% Assays(object)) {
    DefaultAssay(object) <- "RNA"
  }
  
  .save_outputs(
    object = object,
    outdir = outdir,
    method = method,
    batch_col = batch_col,
    dims_use = dims_use,
    resolution = resolution,
    nfeatures = nfeatures,
    npcs = npcs,
    seed = seed,
    n_objects = length(objects)
  )
}

# ============================================================
# Optional Snakemake entry point
# ============================================================

if (exists("snakemake", inherits = FALSE)) {
  run_integration(
    rds_paths = snakemake@input[["rds_list"]],
    method = snakemake@params[["method"]] %||% "Harmony",
    batch_col = snakemake@params[["batch_col"]] %||% "status",
    dims = snakemake@params[["dims"]] %||% 30,
    resolution = snakemake@params[["resolution"]] %||% 0.5,
    nfeatures = snakemake@params[["nfeatures"]] %||% 3000,
    npcs = snakemake@params[["npcs"]] %||% 50,
    seed = snakemake@params[["seed"]] %||% 1234,
    outdir = dirname(snakemake@output[["rds"]])
  )
}