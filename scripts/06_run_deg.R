#!/usr/bin/env Rscript

# ============================================================
# scripts/06_run_deg.R
#
# Reusable differential-expression module for the 10x pipeline.
#
# Supported methods:
#   - pseudo_bulk
#   - FindMarkers
#   - both
#
# The file defines run_deg() and helper functions only.
# No analysis is executed when this file is sourced.
# ============================================================


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})


# ============================================================
# General helpers
# ============================================================

.deg_stop <- function(...) {
  stop(
    paste(...),
    call. = FALSE
  )
}


.deg_message <- function(...) {
  cat(
    "[INFO]",
    ...,
    "\n"
  )
}


.deg_warn <- function(...) {
  cat(
    "[WARN]",
    ...,
    "\n"
  )
}


.deg_safe_name <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  if (!nzchar(x)) {
    x <- "contrast"
  }
  
  x
}


.deg_value <- function(
    row,
    column,
    default = NA_character_
) {
  
  if (
    !column %in% colnames(row) ||
    nrow(row) == 0L
  ) {
    return(default)
  }
  
  value <- row[[column]][[1L]]
  
  if (
    is.null(value) ||
    length(value) == 0L ||
    is.na(value)
  ) {
    return(default)
  }
  
  value <- trimws(
    as.character(value)
  )
  
  if (!nzchar(value)) {
    return(default)
  }
  
  value
}


.deg_validate_numeric <- function(
    value,
    name,
    lower = -Inf,
    upper = Inf
) {
  
  value <- as.numeric(value)
  
  if (
    length(value) != 1L ||
    is.na(value) ||
    value < lower ||
    value > upper
  ) {
    .deg_stop(
      name,
      "must be a numeric value between",
      lower,
      "and",
      upper
    )
  }
  
  value
}


.deg_get_assay_data <- function(
    object,
    assay = "RNA",
    layer = "counts"
) {
  
  if (!assay %in% names(object@assays)) {
    .deg_stop(
      "Assay not found:",
      assay
    )
  }
  
  Seurat::DefaultAssay(object) <- assay
  
  tryCatch(
    Seurat::GetAssayData(
      object,
      assay = assay,
      layer = layer
    ),
    error = function(e) {
      
      tryCatch(
        Seurat::GetAssayData(
          object,
          assay = assay,
          slot = layer
        ),
        error = function(e2) {
          .deg_stop(
            "Unable to read assay",
            assay,
            "layer/slot",
            layer,
            ":",
            conditionMessage(e2)
          )
        }
      )
    }
  )
}


.deg_has_normalized_data <- function(
    object,
    assay = "RNA"
) {
  
  data_mat <- tryCatch(
    .deg_get_assay_data(
      object,
      assay = assay,
      layer = "data"
    ),
    error = function(e) NULL
  )
  
  !is.null(data_mat) &&
    nrow(data_mat) > 0L &&
    ncol(data_mat) > 0L &&
    length(data_mat@x) > 0L
}


.deg_prepare_normalized_object <- function(
    object,
    assay = "RNA"
) {
  
  Seurat::DefaultAssay(object) <- assay
  
  if (
    !.deg_has_normalized_data(
      object,
      assay = assay
    )
  ) {
    
    .deg_message(
      "Normalized data not found;",
      "running NormalizeData()."
    )
    
    object <- Seurat::NormalizeData(
      object,
      assay = assay,
      verbose = FALSE
    )
  }
  
  object
}


# ============================================================
# Contrast validation
# ============================================================

.deg_validate_contrasts <- function(contrasts) {
  
  required_columns <- c(
    "contrast_name",
    "group_col",
    "group1",
    "group2"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(contrasts)
  )
  
  if (length(missing_columns) > 0L) {
    .deg_stop(
      "Contrasts file is missing required column(s):",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  if (nrow(contrasts) == 0L) {
    .deg_stop(
      "The contrasts file contains no contrasts."
    )
  }
  
  for (column in required_columns) {
    
    values <- trimws(
      as.character(
        contrasts[[column]]
      )
    )
    
    if (
      anyNA(values) ||
      any(!nzchar(values))
    ) {
      .deg_stop(
        "Missing or blank values detected in contrasts column:",
        column
      )
    }
  }
  
  duplicated_names <- unique(
    contrasts$contrast_name[
      duplicated(
        contrasts$contrast_name
      )
    ]
  )
  
  if (length(duplicated_names) > 0L) {
    .deg_stop(
      "Duplicated contrast_name values:",
      paste(
        duplicated_names,
        collapse = ", "
      )
    )
  }
  
  optional_columns <- c(
    "subset_col",
    "subset_val",
    "cell_type_col"
  )
  
  for (column in optional_columns) {
    if (!column %in% colnames(contrasts)) {
      contrasts[[column]] <- NA_character_
    }
  }
  
  contrasts
}


.deg_validate_contrast_metadata <- function(
    object,
    group_col,
    group1,
    group2
) {
  
  metadata <- object@meta.data
  
  if (!group_col %in% colnames(metadata)) {
    .deg_stop(
      "Grouping column not found in Seurat metadata:",
      group_col
    )
  }
  
  observed_groups <- unique(
    as.character(
      metadata[[group_col]]
    )
  )
  
  missing_groups <- setdiff(
    c(group1, group2),
    observed_groups
  )
  
  if (length(missing_groups) > 0L) {
    .deg_stop(
      "Group level(s) not found in",
      group_col,
      ":",
      paste(
        missing_groups,
        collapse = ", "
      )
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# Optional cell subsetting
# ============================================================

.deg_subset_object <- function(
    object,
    subset_col = NA_character_,
    subset_val = NA_character_
) {
  
  has_subset_col <- !is.na(subset_col) &&
    nzchar(subset_col)
  
  has_subset_val <- !is.na(subset_val) &&
    nzchar(subset_val)
  
  if (!has_subset_col && !has_subset_val) {
    return(object)
  }
  
  if (xor(has_subset_col, has_subset_val)) {
    .deg_stop(
      "Both subset_col and subset_val must be provided."
    )
  }
  
  if (!subset_col %in% colnames(object@meta.data)) {
    .deg_stop(
      "subset_col not found in Seurat metadata:",
      subset_col
    )
  }
  
  keep <- !is.na(
    object@meta.data[[subset_col]]
  ) &
    as.character(
      object@meta.data[[subset_col]]
    ) == subset_val
  
  selected_cells <- rownames(
    object@meta.data
  )[keep]
  
  if (length(selected_cells) == 0L) {
    .deg_stop(
      "No cells remain after subsetting",
      subset_col,
      "=",
      subset_val
    )
  }
  
  object <- subset(
    object,
    cells = selected_cells
  )
  
  .deg_message(
    "Subset:",
    subset_col,
    "=",
    subset_val,
    "->",
    ncol(object),
    "cells"
  )
  
  object
}


# ============================================================
# Pseudo-bulk aggregation
# ============================================================

.deg_aggregate_pseudobulk <- function(
    object,
    donor_col,
    group_col,
    group1,
    group2,
    assay = "RNA"
) {
  
  metadata <- object@meta.data
  
  required_metadata <- c(
    donor_col,
    group_col
  )
  
  missing_metadata <- setdiff(
    required_metadata,
    colnames(metadata)
  )
  
  if (length(missing_metadata) > 0L) {
    .deg_stop(
      "Pseudo-bulk metadata column(s) not found:",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  }
  
  keep <- !is.na(
    metadata[[donor_col]]
  ) &
    !is.na(
      metadata[[group_col]]
    ) &
    as.character(
      metadata[[group_col]]
    ) %in% c(group1, group2)
  
  metadata <- metadata[
    keep,
    ,
    drop = FALSE
  ]
  
  if (nrow(metadata) == 0L) {
    .deg_stop(
      "No cells remain for pseudo-bulk aggregation."
    )
  }
  
  counts <- .deg_get_assay_data(
    object,
    assay = assay,
    layer = "counts"
  )
  
  counts <- counts[
    ,
    rownames(metadata),
    drop = FALSE
  ]
  
  metadata$.deg_donor <- as.character(
    metadata[[donor_col]]
  )
  
  metadata$.deg_group <- as.character(
    metadata[[group_col]]
  )
  
  metadata$.deg_pb_id <- paste(
    metadata$.deg_donor,
    metadata$.deg_group,
    sep = "__"
  )
  
  pb_ids <- unique(
    metadata$.deg_pb_id
  )
  
  aggregated <- lapply(
    pb_ids,
    function(pb_id) {
      
      cells <- rownames(metadata)[
        metadata$.deg_pb_id == pb_id
      ]
      
      Matrix::rowSums(
        counts[
          ,
          cells,
          drop = FALSE
        ]
      )
    }
  )
  
  pb_counts <- do.call(
    cbind,
    aggregated
  )
  
  rownames(pb_counts) <- rownames(counts)
  colnames(pb_counts) <- pb_ids
  
  split_ids <- strsplit(
    pb_ids,
    "__",
    fixed = TRUE
  )
  
  pb_meta <- data.frame(
    donor = vapply(
      split_ids,
      `[[`,
      character(1L),
      1L
    ),
    group = vapply(
      split_ids,
      `[[`,
      character(1L),
      2L
    ),
    row.names = pb_ids,
    stringsAsFactors = FALSE
  )
  
  colnames(pb_meta) <- c(
    donor_col,
    group_col
  )
  
  list(
    counts = pb_counts,
    metadata = pb_meta
  )
}


# ============================================================
# Pseudo-bulk DESeq2
# ============================================================

.deg_run_pseudo_bulk <- function(
    object,
    donor_col,
    group_col,
    group1,
    group2,
    fdr_thr,
    assay = "RNA",
    min_replicates = 2L,
    min_total_count = 10L
) {
  
  if (!requireNamespace(
    "DESeq2",
    quietly = TRUE
  )) {
    .deg_stop(
      "Package 'DESeq2' is required",
      "for pseudo-bulk analysis."
    )
  }
  
  if (!requireNamespace(
    "Matrix",
    quietly = TRUE
  )) {
    .deg_stop(
      "Package 'Matrix' is required",
      "for pseudo-bulk analysis."
    )
  }
  
  aggregated <- .deg_aggregate_pseudobulk(
    object = object,
    donor_col = donor_col,
    group_col = group_col,
    group1 = group1,
    group2 = group2,
    assay = assay
  )
  
  pb_counts <- aggregated$counts
  pb_meta <- aggregated$metadata
  
  replicate_counts <- table(
    pb_meta[[group_col]]
  )
  
  missing_groups <- setdiff(
    c(group1, group2),
    names(replicate_counts)
  )
  
  if (length(missing_groups) > 0L) {
    .deg_stop(
      "Pseudo-bulk group(s) missing:",
      paste(
        missing_groups,
        collapse = ", "
      )
    )
  }
  
  insufficient_groups <- names(
    replicate_counts
  )[
    replicate_counts < min_replicates
  ]
  
  if (length(insufficient_groups) > 0L) {
    .deg_stop(
      "Pseudo-bulk requires at least",
      min_replicates,
      "replicates per group. Insufficient group(s):",
      paste(
        insufficient_groups,
        collapse = ", "
      )
    )
  }
  
  keep_genes <- Matrix::rowSums(
    pb_counts
  ) >= min_total_count
  
  pb_counts <- pb_counts[
    keep_genes,
    ,
    drop = FALSE
  ]
  
  if (nrow(pb_counts) == 0L) {
    .deg_stop(
      "No genes remain after pseudo-bulk count filtering."
    )
  }
  
  pb_meta[[group_col]] <- factor(
    pb_meta[[group_col]],
    levels = c(
      group2,
      group1
    )
  )
  
  design_formula <- stats::as.formula(
    paste0(
      "~ `",
      group_col,
      "`"
    )
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = pb_meta,
    design = design_formula
  )
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  result <- DESeq2::results(
    dds,
    contrast = c(
      group_col,
      group1,
      group2
    ),
    alpha = fdr_thr
  )
  
  result <- as.data.frame(
    result
  )
  
  result$gene <- rownames(result)
  
  result <- result |>
    dplyr::transmute(
      gene = gene,
      baseMean = baseMean,
      avg_log2FC = log2FoldChange,
      lfcSE = lfcSE,
      stat = stat,
      p_val = pvalue,
      p_val_adj = padj
    ) |>
    dplyr::arrange(
      is.na(p_val_adj),
      p_val_adj,
      dplyr::desc(
        abs(avg_log2FC)
      )
    )
  
  attr(
    result,
    "dds"
  ) <- dds
  
  result
}


# ============================================================
# Seurat FindMarkers
# ============================================================

.deg_run_find_markers <- function(
    object,
    group_col,
    group1,
    group2,
    lfc_thr,
    assay = "RNA",
    test_use = "wilcox",
    min_pct = 0.1
) {
  
  object <- .deg_prepare_normalized_object(
    object,
    assay = assay
  )
  
  object@meta.data[[group_col]] <- as.character(
    object@meta.data[[group_col]]
  )
  
  Seurat::Idents(object) <- group_col
  
  result <- Seurat::FindMarkers(
    object = object,
    ident.1 = group1,
    ident.2 = group2,
    assay = assay,
    test.use = test_use,
    logfc.threshold = lfc_thr,
    min.pct = min_pct,
    verbose = FALSE
  )
  
  if (nrow(result) == 0L) {
    return(
      data.frame()
    )
  }
  
  result$gene <- rownames(result)
  
  if (
    !"avg_log2FC" %in% colnames(result) &&
    "avg_logFC" %in% colnames(result)
  ) {
    result$avg_log2FC <- result$avg_logFC
  }
  
  required_columns <- c(
    "gene",
    "avg_log2FC",
    "p_val",
    "p_val_adj"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(result)
  )
  
  if (length(missing_columns) > 0L) {
    .deg_stop(
      "FindMarkers result is missing column(s):",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  result |>
    dplyr::relocate(
      gene
    ) |>
    dplyr::arrange(
      is.na(p_val_adj),
      p_val_adj,
      dplyr::desc(
        abs(avg_log2FC)
      )
    )
}


# ============================================================
# Volcano plot
# ============================================================

.deg_make_volcano <- function(
    deg,
    contrast_name,
    group1,
    group2,
    lfc_thr,
    fdr_thr,
    top_n
) {
  
  plot_data <- deg |>
    dplyr::filter(
      !is.na(avg_log2FC),
      !is.na(p_val_adj)
    ) |>
    dplyr::mutate(
      regulation = dplyr::case_when(
        avg_log2FC >= lfc_thr &
          p_val_adj < fdr_thr ~ "Up",
        
        avg_log2FC <= -lfc_thr &
          p_val_adj < fdr_thr ~ "Down",
        
        TRUE ~ "NS"
      ),
      
      neg_log10_padj = -log10(
        pmax(
          p_val_adj,
          .Machine$double.xmin
        )
      )
    )
  
  if (nrow(plot_data) == 0L) {
    return(NULL)
  }
  
  ranked_labels <- plot_data |>
    dplyr::filter(
      regulation != "NS"
    ) |>
    dplyr::arrange(
      p_val_adj,
      dplyr::desc(
        abs(avg_log2FC)
      )
    ) |>
    dplyr::slice_head(
      n = top_n
    ) |>
    dplyr::pull(
      gene
    )
  
  plot_data$label <- ifelse(
    plot_data$gene %in% ranked_labels,
    plot_data$gene,
    NA_character_
  )
  
  counts <- table(
    factor(
      plot_data$regulation,
      levels = c(
        "Up",
        "Down",
        "NS"
      )
    )
  )
  
  palette <- c(
    Up = "#740001",
    Down = "#6497B1",
    NS = "grey70"
  )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = avg_log2FC,
      y = neg_log10_padj,
      colour = regulation
    )
  ) +
    ggplot2::geom_point(
      size = 0.9,
      alpha = 0.75
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        label = label
      ),
      size = 2.7,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.2,
      colour = "grey20",
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -lfc_thr,
        lfc_thr
      ),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(
        fdr_thr
      ),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    ggplot2::scale_colour_manual(
      values = palette,
      breaks = c(
        "Up",
        "Down",
        "NS"
      )
    ) +
    ggplot2::labs(
      title = contrast_name,
      subtitle = paste0(
        group1,
        " vs ",
        group2,
        " | Up: ",
        counts[["Up"]],
        " | Down: ",
        counts[["Down"]]
      ),
      x = expression(
        Log[2] ~ fold ~ change
      ),
      y = expression(
        -Log[10] ~ adjusted ~ italic(P)
      ),
      colour = "Regulation"
    ) +
    ggplot2::theme_classic(
      base_size = 10
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9,
        colour = "grey35"
      ),
      legend.position = "right"
    )
}


# ============================================================
# Heatmap
# ============================================================

.deg_make_heatmap <- function(
    object,
    deg,
    group_col,
    group1,
    group2,
    contrast_name,
    fdr_thr,
    lfc_thr,
    top_n,
    assay = "RNA"
) {
  
  if (!requireNamespace(
    "ComplexHeatmap",
    quietly = TRUE
  )) {
    .deg_warn(
      "Package 'ComplexHeatmap' is unavailable;",
      "heatmap skipped."
    )
    
    return(NULL)
  }
  
  if (!requireNamespace(
    "circlize",
    quietly = TRUE
  )) {
    .deg_warn(
      "Package 'circlize' is unavailable;",
      "heatmap skipped."
    )
    
    return(NULL)
  }
  
  significant <- deg |>
    dplyr::filter(
      !is.na(p_val_adj),
      !is.na(avg_log2FC),
      p_val_adj < fdr_thr,
      abs(avg_log2FC) >= lfc_thr
    ) |>
    dplyr::arrange(
      p_val_adj,
      dplyr::desc(
        abs(avg_log2FC)
      )
    ) |>
    dplyr::slice_head(
      n = top_n
    )
  
  top_genes <- unique(
    significant$gene
  )
  
  if (length(top_genes) < 3L) {
    .deg_warn(
      "Too few significant genes for heatmap:",
      contrast_name
    )
    
    return(NULL)
  }
  
  object <- .deg_prepare_normalized_object(
    object,
    assay = assay
  )
  
  data_mat <- .deg_get_assay_data(
    object,
    assay = assay,
    layer = "data"
  )
  
  top_genes <- intersect(
    top_genes,
    rownames(data_mat)
  )
  
  if (length(top_genes) < 3L) {
    .deg_warn(
      "Too few selected genes are present in assay",
      assay,
      "for:",
      contrast_name
    )
    
    return(NULL)
  }
  
  metadata <- object@meta.data
  
  keep_cells <- !is.na(
    metadata[[group_col]]
  ) &
    as.character(
      metadata[[group_col]]
    ) %in% c(
      group1,
      group2
    )
  
  cells <- rownames(metadata)[
    keep_cells
  ]
  
  metadata <- metadata[
    cells,
    ,
    drop = FALSE
  ]
  
  group_factor <- factor(
    as.character(
      metadata[[group_col]]
    ),
    levels = c(
      group2,
      group1
    )
  )
  
  column_order <- order(
    group_factor
  )
  
  cells <- cells[
    column_order
  ]
  
  group_factor <- group_factor[
    column_order
  ]
  
  heatmap_matrix <- as.matrix(
    data_mat[
      top_genes,
      cells,
      drop = FALSE
    ]
  )
  
  row_sd <- apply(
    heatmap_matrix,
    1L,
    stats::sd,
    na.rm = TRUE
  )
  
  heatmap_matrix <- heatmap_matrix[
    !is.na(row_sd) &
      row_sd > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(heatmap_matrix) < 3L) {
    .deg_warn(
      "Too few variable genes for heatmap:",
      contrast_name
    )
    
    return(NULL)
  }
  
  heatmap_matrix <- t(
    scale(
      t(
        heatmap_matrix
      )
    )
  )
  
  heatmap_matrix[
    !is.finite(
      heatmap_matrix
    )
  ] <- 0
  
  heatmap_matrix <- pmax(
    pmin(
      heatmap_matrix,
      2
    ),
    -2
  )
  
  annotation_colors <- c(
    group2 = "#4DBBD5",
    group1 = "#D55E00"
  )
  
  names(annotation_colors) <- c(
    group2,
    group1
  )
  
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = group_factor,
    col = list(
      Group = annotation_colors
    ),
    annotation_name_side = "left"
  )
  
  ComplexHeatmap::Heatmap(
    heatmap_matrix,
    name = "Z-score",
    col = circlize::colorRamp2(
      c(
        -2,
        0,
        2
      ),
      c(
        "#6497B1",
        "#F7F7F7",
        "#740001"
      )
    ),
    top_annotation = top_annotation,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(
      fontsize = 7
    ),
    column_title = contrast_name,
    column_title_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
    use_raster = ncol(
      heatmap_matrix
    ) > 1000
  )
}


# ============================================================
# Pathway enrichment
# ============================================================

.deg_get_orgdb <- function(species) {
  
  species <- tolower(
    trimws(
      species
    )
  )
  
  if (
    species %in% c(
      "human",
      "homo sapiens",
      "hs",
      "hsa"
    )
  ) {
    
    if (!requireNamespace(
      "org.Hs.eg.db",
      quietly = TRUE
    )) {
      .deg_stop(
        "Package 'org.Hs.eg.db' is required",
        "for human pathway enrichment."
      )
    }
    
    return(
      org.Hs.eg.db::org.Hs.eg.db
    )
  }
  
  if (
    species %in% c(
      "mouse",
      "mus musculus",
      "mm",
      "mmu"
    )
  ) {
    
    if (!requireNamespace(
      "org.Mm.eg.db",
      quietly = TRUE
    )) {
      .deg_stop(
        "Package 'org.Mm.eg.db' is required",
        "for mouse pathway enrichment."
      )
    }
    
    return(
      org.Mm.eg.db::org.Mm.eg.db
    )
  }
  
  .deg_stop(
    "Unsupported species for pathway enrichment:",
    species
  )
}


.deg_run_pathway <- function(
    deg,
    contrast_name,
    species,
    fdr_thr,
    lfc_thr,
    direction = "Up",
    min_genes = 10L
) {
  
  if (!requireNamespace(
    "clusterProfiler",
    quietly = TRUE
  )) {
    .deg_warn(
      "Package 'clusterProfiler' is unavailable;",
      "pathway enrichment skipped."
    )
    
    return(NULL)
  }
  
  direction <- match.arg(
    direction,
    choices = c(
      "Up",
      "Down"
    )
  )
  
  selected <- if (direction == "Up") {
    
    deg |>
      dplyr::filter(
        !is.na(p_val_adj),
        !is.na(avg_log2FC),
        p_val_adj < fdr_thr,
        avg_log2FC >= lfc_thr
      )
    
  } else {
    
    deg |>
      dplyr::filter(
        !is.na(p_val_adj),
        !is.na(avg_log2FC),
        p_val_adj < fdr_thr,
        avg_log2FC <= -lfc_thr
      )
  }
  
  genes <- unique(
    selected$gene
  )
  
  if (length(genes) < min_genes) {
    .deg_warn(
      "Too few",
      tolower(direction),
      "genes for pathway enrichment in:",
      contrast_name
    )
    
    return(NULL)
  }
  
  org_db <- .deg_get_orgdb(
    species
  )
  
  mapped_ids <- tryCatch(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org_db
    ),
    error = function(e) {
      .deg_warn(
        "Gene identifier conversion failed:",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (
    is.null(mapped_ids) ||
    nrow(mapped_ids) == 0L
  ) {
    return(NULL)
  }
  
  enrichment <- clusterProfiler::enrichGO(
    gene = unique(
      mapped_ids$ENTREZID
    ),
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  enrichment_table <- as.data.frame(
    enrichment
  )
  
  if (nrow(enrichment_table) == 0L) {
    return(NULL)
  }
  
  enrichment_plot <- clusterProfiler::dotplot(
    enrichment,
    showCategory = 15,
    title = paste0(
      "GO:BP — ",
      contrast_name,
      " — ",
      direction
    )
  ) +
    ggplot2::theme_classic(
      base_size = 9
    )
  
  list(
    result = enrichment,
    table = enrichment_table,
    plot = enrichment_plot,
    direction = direction
  )
}


# ============================================================
# HTML report
# ============================================================

.deg_generate_report <- function(
    results,
    report_file,
    tables_dir,
    plots_dir
) {
  
  if (!requireNamespace(
    "rmarkdown",
    quietly = TRUE
  )) {
    .deg_warn(
      "Package 'rmarkdown' is unavailable;",
      "HTML report skipped."
    )
    
    return(NULL)
  }
  
  if (!requireNamespace(
    "knitr",
    quietly = TRUE
  )) {
    .deg_warn(
      "Package 'knitr' is unavailable;",
      "HTML report skipped."
    )
    
    return(NULL)
  }
  
  report_dir <- dirname(
    report_file
  )
  
  dir.create(
    report_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  report_rmd <- file.path(
    report_dir,
    "deg_report.Rmd"
  )
  
  lines <- c(
    "---",
    "title: \"Differential Expression Analysis\"",
    "date: \"`r format(Sys.time(), '%Y-%m-%d %H:%M')`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: flatly",
    "    code_folding: hide",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(",
    "  echo = FALSE,",
    "  warning = FALSE,",
    "  message = FALSE",
    ")",
    "```",
    "",
    "# Summary",
    "",
    paste0(
      "Number of completed contrasts: **",
      length(results),
      "**"
    ),
    ""
  )
  
  for (contrast_name in names(results)) {
    
    safe_name <- .deg_safe_name(
      contrast_name
    )
    
    table_file <- normalizePath(
      file.path(
        tables_dir,
        paste0(
          safe_name,
          "_DEG.tsv"
        )
      ),
      mustWork = FALSE
    )
    
    volcano_file <- normalizePath(
      file.path(
        plots_dir,
        paste0(
          safe_name,
          "_volcano.png"
        )
      ),
      mustWork = FALSE
    )
    
    heatmap_file <- normalizePath(
      file.path(
        plots_dir,
        paste0(
          safe_name,
          "_heatmap.png"
        )
      ),
      mustWork = FALSE
    )
    
    lines <- c(
      lines,
      paste0(
        "# ",
        contrast_name
      ),
      "",
      "## DEG table",
      "",
      "```{r}",
      paste0(
        "deg_table <- read.delim(",
        deparse(table_file),
        ", check.names = FALSE)"
      ),
      "knitr::kable(",
      "  head(deg_table, 50),",
      "  digits = 4,",
      "  caption = \"Top differential-expression results\"",
      ")",
      "```",
      ""
    )
    
    if (file.exists(volcano_file)) {
      lines <- c(
        lines,
        "## Volcano plot",
        "",
        "```{r}",
        paste0(
          "knitr::include_graphics(",
          deparse(volcano_file),
          ")"
        ),
        "```",
        ""
      )
    }
    
    if (file.exists(heatmap_file)) {
      lines <- c(
        lines,
        "## Heatmap",
        "",
        "```{r}",
        paste0(
          "knitr::include_graphics(",
          deparse(heatmap_file),
          ")"
        ),
        "```",
        ""
      )
    }
  }
  
  writeLines(
    lines,
    report_rmd
  )
  
  rendered <- tryCatch(
    rmarkdown::render(
      input = report_rmd,
      output_file = basename(
        report_file
      ),
      output_dir = report_dir,
      quiet = TRUE,
      envir = new.env(
        parent = globalenv()
      )
    ),
    error = function(e) {
      .deg_warn(
        "HTML report generation failed:",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  rendered
}


# ============================================================
# Main reusable function
# ============================================================

run_deg <- function(
    rds_in,
    contrasts_file,
    method = "pseudo_bulk",
    donor_col = "donor_id",
    fdr_thr = 0.05,
    lfc_thr = 0.5,
    top_n = 20,
    run_pathway = TRUE,
    species = "human",
    outdir = "results/06_deg",
    assay = "RNA",
    findmarkers_test = "wilcox",
    min_pct = 0.1,
    min_pseudobulk_replicates = 2L,
    min_pseudobulk_count = 10L,
    generate_report = TRUE,
    seed = 1234
) {
  
  # ----------------------------------------------------------
  # Validate parameters
  # ----------------------------------------------------------
  
  if (
    is.null(rds_in) ||
    length(rds_in) != 1L ||
    is.na(rds_in) ||
    !file.exists(rds_in)
  ) {
    .deg_stop(
      "Integrated RDS not found:",
      rds_in
    )
  }
  
  if (
    is.null(contrasts_file) ||
    length(contrasts_file) != 1L ||
    is.na(contrasts_file) ||
    !file.exists(contrasts_file)
  ) {
    .deg_stop(
      "Contrasts file not found:",
      contrasts_file
    )
  }
  
  method <- tolower(
    trimws(
      method
    )
  )
  
  method_aliases <- c(
    pseudobulk = "pseudo_bulk",
    pseudo_bulk = "pseudo_bulk",
    findmarkers = "findmarkers",
    find_markers = "findmarkers",
    both = "both"
  )
  
  if (!method %in% names(method_aliases)) {
    .deg_stop(
      "Unsupported DEG method:",
      method,
      ". Use pseudo_bulk, FindMarkers or both."
    )
  }
  
  method <- unname(
    method_aliases[[method]]
  )
  
  fdr_thr <- .deg_validate_numeric(
    fdr_thr,
    "fdr_thr",
    lower = 0,
    upper = 1
  )
  
  lfc_thr <- .deg_validate_numeric(
    lfc_thr,
    "lfc_thr",
    lower = 0
  )
  
  top_n <- as.integer(
    top_n
  )
  
  if (
    is.na(top_n) ||
    top_n < 1L
  ) {
    .deg_stop(
      "top_n must be a positive integer."
    )
  }
  
  min_pseudobulk_replicates <- as.integer(
    min_pseudobulk_replicates
  )
  
  min_pseudobulk_count <- as.integer(
    min_pseudobulk_count
  )
  
  set.seed(
    as.integer(seed)
  )
  
  # ----------------------------------------------------------
  # Prepare output directories
  # ----------------------------------------------------------
  
  tables_dir <- file.path(
    outdir,
    "tables"
  )
  
  plots_dir <- file.path(
    outdir,
    "plots"
  )
  
  pathways_dir <- file.path(
    outdir,
    "pathways"
  )
  
  objects_dir <- file.path(
    outdir,
    "objects"
  )
  
  dir.create(
    tables_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    plots_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    objects_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (isTRUE(run_pathway)) {
    dir.create(
      pathways_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  # ----------------------------------------------------------
  # Read inputs
  # ----------------------------------------------------------
  
  .deg_message(
    "Reading integrated object:",
    rds_in
  )
  
  object <- readRDS(
    rds_in
  )
  
  if (!inherits(
    object,
    "Seurat"
  )) {
    .deg_stop(
      "Integrated RDS does not contain a Seurat object."
    )
  }
  
  if (!assay %in% names(object@assays)) {
    .deg_stop(
      "Requested assay is absent from the Seurat object:",
      assay
    )
  }
  
  .deg_message(
    "Reading contrasts:",
    contrasts_file
  )
  
  contrasts <- read.csv(
    contrasts_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  contrasts <- .deg_validate_contrasts(
    contrasts
  )
  
  .deg_message(
    "DEG method:",
    method
  )
  
  .deg_message(
    "Donor column:",
    donor_col
  )
  
  .deg_message(
    "FDR threshold:",
    fdr_thr
  )
  
  .deg_message(
    "log2FC threshold:",
    lfc_thr
  )
  
  .deg_message(
    "Contrasts:",
    nrow(contrasts)
  )
  
  # ----------------------------------------------------------
  # Run contrasts
  # ----------------------------------------------------------
  
  all_results <- list()
  status_rows <- list()
  
  for (i in seq_len(
    nrow(contrasts)
  )) {
    
    contrast_row <- contrasts[
      i,
      ,
      drop = FALSE
    ]
    
    contrast_name <- .deg_value(
      contrast_row,
      "contrast_name"
    )
    
    group_col <- .deg_value(
      contrast_row,
      "group_col"
    )
    
    group1 <- .deg_value(
      contrast_row,
      "group1"
    )
    
    group2 <- .deg_value(
      contrast_row,
      "group2"
    )
    
    subset_col <- .deg_value(
      contrast_row,
      "subset_col",
      NA_character_
    )
    
    subset_val <- .deg_value(
      contrast_row,
      "subset_val",
      NA_character_
    )
    
    safe_name <- .deg_safe_name(
      contrast_name
    )
    
    cat(
      "\n[INFO] ─────────────────────────────────────────────\n"
    )
    
    .deg_message(
      "Contrast:",
      contrast_name
    )
    
    .deg_message(
      "Comparison:",
      group1,
      "vs",
      group2
    )
    
    contrast_result <- tryCatch({
      
      object_subset <- .deg_subset_object(
        object = object,
        subset_col = subset_col,
        subset_val = subset_val
      )
      
      .deg_validate_contrast_metadata(
        object = object_subset,
        group_col = group_col,
        group1 = group1,
        group2 = group2
      )
      
      deg <- NULL
      method_used <- NULL
      pseudobulk_error <- NULL
      
      if (
        method %in% c(
          "pseudo_bulk",
          "both"
        )
      ) {
        
        .deg_message(
          "Running pseudo-bulk DESeq2."
        )
        
        deg <- tryCatch(
          .deg_run_pseudo_bulk(
            object = object_subset,
            donor_col = donor_col,
            group_col = group_col,
            group1 = group1,
            group2 = group2,
            fdr_thr = fdr_thr,
            assay = assay,
            min_replicates =
              min_pseudobulk_replicates,
            min_total_count =
              min_pseudobulk_count
          ),
          error = function(e) {
            pseudobulk_error <<- conditionMessage(
              e
            )
            
            NULL
          }
        )
        
        if (!is.null(deg)) {
          method_used <- "pseudo_bulk"
        }
      }
      
      if (
        method == "findmarkers" ||
        is.null(deg)
      ) {
        
        if (!is.null(pseudobulk_error)) {
          .deg_warn(
            "Pseudo-bulk failed:",
            pseudobulk_error
          )
          
          .deg_warn(
            "Falling back to Seurat FindMarkers."
          )
        }
        
        .deg_message(
          "Running Seurat FindMarkers."
        )
        
        deg <- .deg_run_find_markers(
          object = object_subset,
          group_col = group_col,
          group1 = group1,
          group2 = group2,
          lfc_thr = lfc_thr,
          assay = assay,
          test_use = findmarkers_test,
          min_pct = min_pct
        )
        
        method_used <- "FindMarkers"
      }
      
      if (
        is.null(deg) ||
        nrow(deg) == 0L
      ) {
        .deg_stop(
          "No differential-expression results were produced."
        )
      }
      
      deg <- deg |>
        dplyr::mutate(
          contrast_name = contrast_name,
          group_col = group_col,
          group1 = group1,
          group2 = group2,
          method = method_used,
          significant = !is.na(
            p_val_adj
          ) &
            p_val_adj < fdr_thr &
            abs(
              avg_log2FC
            ) >= lfc_thr,
          regulation = dplyr::case_when(
            significant &
              avg_log2FC > 0 ~ "Up",
            
            significant &
              avg_log2FC < 0 ~ "Down",
            
            TRUE ~ "NS"
          )
        ) |>
        dplyr::relocate(
          contrast_name,
          group_col,
          group1,
          group2,
          method,
          gene
        )
      
      table_file <- file.path(
        tables_dir,
        paste0(
          safe_name,
          "_DEG.tsv"
        )
      )
      
      utils::write.table(
        deg,
        table_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        na = ""
      )
      
      .deg_message(
        "Saved DEG table:",
        table_file
      )
      
      # ------------------------------------------------------
      # Volcano plot
      # ------------------------------------------------------
      
      volcano <- .deg_make_volcano(
        deg = deg,
        contrast_name = contrast_name,
        group1 = group1,
        group2 = group2,
        lfc_thr = lfc_thr,
        fdr_thr = fdr_thr,
        top_n = top_n
      )
      
      if (!is.null(volcano)) {
        
        volcano_pdf <- file.path(
          plots_dir,
          paste0(
            safe_name,
            "_volcano.pdf"
          )
        )
        
        volcano_png <- file.path(
          plots_dir,
          paste0(
            safe_name,
            "_volcano.png"
          )
        )
        
        ggplot2::ggsave(
          filename = volcano_pdf,
          plot = volcano,
          width = 7,
          height = 6,
          units = "in",
          device = cairo_pdf
        )
        
        ggplot2::ggsave(
          filename = volcano_png,
          plot = volcano,
          width = 7,
          height = 6,
          units = "in",
          dpi = 300
        )
      }
      
      # ------------------------------------------------------
      # Heatmap
      # ------------------------------------------------------
      
      heatmap <- .deg_make_heatmap(
        object = object_subset,
        deg = deg,
        group_col = group_col,
        group1 = group1,
        group2 = group2,
        contrast_name = contrast_name,
        fdr_thr = fdr_thr,
        lfc_thr = lfc_thr,
        top_n = top_n,
        assay = assay
      )
      
      if (!is.null(heatmap)) {
        
        heatmap_pdf <- file.path(
          plots_dir,
          paste0(
            safe_name,
            "_heatmap.pdf"
          )
        )
        
        heatmap_png <- file.path(
          plots_dir,
          paste0(
            safe_name,
            "_heatmap.png"
          )
        )
        
        grDevices::pdf(
          heatmap_pdf,
          width = 8,
          height = 10
        )
        
        ComplexHeatmap::draw(
          heatmap
        )
        
        grDevices::dev.off()
        
        grDevices::png(
          heatmap_png,
          width = 8,
          height = 10,
          units = "in",
          res = 300
        )
        
        ComplexHeatmap::draw(
          heatmap
        )
        
        grDevices::dev.off()
      }
      
      # ------------------------------------------------------
      # Pathway enrichment
      # ------------------------------------------------------
      
      pathway_results <- list()
      
      if (isTRUE(run_pathway)) {
        
        for (direction in c(
          "Up",
          "Down"
        )) {
          
          pathway <- tryCatch(
            .deg_run_pathway(
              deg = deg,
              contrast_name = contrast_name,
              species = species,
              fdr_thr = fdr_thr,
              lfc_thr = lfc_thr,
              direction = direction
            ),
            error = function(e) {
              .deg_warn(
                "Pathway analysis failed for",
                contrast_name,
                direction,
                ":",
                conditionMessage(e)
              )
              
              NULL
            }
          )
          
          if (!is.null(pathway)) {
            
            direction_safe <- tolower(
              direction
            )
            
            pathway_table <- file.path(
              pathways_dir,
              paste0(
                safe_name,
                "_",
                direction_safe,
                "_GO_BP.tsv"
              )
            )
            
            pathway_plot <- file.path(
              plots_dir,
              paste0(
                safe_name,
                "_",
                direction_safe,
                "_GO_BP.pdf"
              )
            )
            
            utils::write.table(
              pathway$table,
              pathway_table,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE,
              na = ""
            )
            
            ggplot2::ggsave(
              filename = pathway_plot,
              plot = pathway$plot,
              width = 9,
              height = 7,
              units = "in",
              device = cairo_pdf
            )
            
            pathway_results[[
              direction
            ]] <- pathway
          }
        }
      }
      
      n_up <- sum(
        deg$regulation == "Up",
        na.rm = TRUE
      )
      
      n_down <- sum(
        deg$regulation == "Down",
        na.rm = TRUE
      )
      
      list(
        deg = deg,
        volcano = volcano,
        heatmap = heatmap,
        pathways = pathway_results,
        method = method_used,
        n_up = n_up,
        n_down = n_down
      )
      
    }, error = function(e) {
      
      .deg_warn(
        "Contrast failed:",
        contrast_name,
        "—",
        conditionMessage(e)
      )
      
      attr(
        list(),
        "error"
      ) <- conditionMessage(e)
      
      NULL
    })
    
    if (is.null(contrast_result)) {
      
      status_rows[[length(
        status_rows
      ) + 1L]] <- data.frame(
        contrast_name = contrast_name,
        status = "FAILED",
        method = NA_character_,
        n_up = NA_integer_,
        n_down = NA_integer_,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    all_results[[
      contrast_name
    ]] <- contrast_result
    
    status_rows[[length(
      status_rows
    ) + 1L]] <- data.frame(
      contrast_name = contrast_name,
      status = "SUCCESS",
      method = contrast_result$method,
      n_up = contrast_result$n_up,
      n_down = contrast_result$n_down,
      stringsAsFactors = FALSE
    )
    
    invisible(
      gc()
    )
  }
  
  # ----------------------------------------------------------
  # Save summary
  # ----------------------------------------------------------
  
  summary_table <- if (
    length(status_rows) > 0L
  ) {
    do.call(
      rbind,
      status_rows
    )
  } else {
    data.frame(
      contrast_name = character(),
      status = character(),
      method = character(),
      n_up = integer(),
      n_down = integer()
    )
  }
  
  summary_file <- file.path(
    outdir,
    "deg_summary.tsv"
  )
  
  utils::write.table(
    summary_table,
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  
  results_file <- file.path(
    objects_dir,
    "deg_results.rds"
  )
  
  saveRDS(
    all_results,
    results_file
  )
  
  # ----------------------------------------------------------
  # Generate report
  # ----------------------------------------------------------
  
  report_file <- file.path(
    outdir,
    "deg_report.html"
  )
  
  rendered_report <- NULL
  
  if (
    isTRUE(generate_report) &&
    length(all_results) > 0L
  ) {
    rendered_report <- .deg_generate_report(
      results = all_results,
      report_file = report_file,
      tables_dir = tables_dir,
      plots_dir = plots_dir
    )
  }
  
  # ----------------------------------------------------------
  # Final status
  # ----------------------------------------------------------
  
  successful <- sum(
    summary_table$status == "SUCCESS"
  )
  
  failed <- sum(
    summary_table$status == "FAILED"
  )
  
  .deg_message(
    "DEG analysis complete."
  )
  
  .deg_message(
    "Successful contrasts:",
    successful
  )
  
  .deg_message(
    "Failed contrasts:",
    failed
  )
  
  .deg_message(
    "Output directory:",
    outdir
  )
  
  invisible(
    list(
      results = all_results,
      summary = summary_table,
      report = rendered_report,
      results_rds = results_file,
      output_dir = outdir,
      tables_dir = tables_dir,
      plots_dir = plots_dir,
      pathways_dir = if (
        isTRUE(run_pathway)
      ) {
        pathways_dir
      } else {
        NULL
      }
    )
  )
}