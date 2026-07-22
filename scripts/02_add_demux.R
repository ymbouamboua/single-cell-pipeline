#!/usr/bin/env Rscript

# ============================================================ #
# scripts/02_add_demux.R
#
# Add Demuxafy and/or bulk-genotype assignments to a Seurat
# object.
#
# Main output:
#   results/02_demuxed/<sample>_demuxed.rds
#
# Standardized metadata:
#   status                 Final donor/sample assignment
#   demux_status           Demuxafy consensus assignment
#   demux_class            Singlet/doublet/unassigned class
#   BulkSample             Bulk-genotype assignment
#   demux_is_doublet       Logical doublet flag
#   demux_is_unassigned    Logical unassigned flag
#   demux_*                Original Demuxafy columns
#
# This file defines run_demux(). It only executes automatically
# when called by Snakemake with a real `snakemake` object.
# ============================================================ #


suppressPackageStartupMessages({
  library(Seurat)
})


# ============================================================ #
# General helpers
# ============================================================ #

`%||%` <- function(x, y) {
  
  if (
    is.null(x) ||
    length(x) == 0L
  ) {
    return(y)
  }
  
  if (
    length(x) == 1L &&
    is.na(x)
  ) {
    return(y)
  }
  
  x
}


.msg <- function(verbose = TRUE) {
  
  function(...) {
    
    if (isTRUE(verbose)) {
      cat(
        "[INFO]",
        ...,
        "\n"
      )
    }
  }
}


.warn <- function(...) {
  
  cat(
    "[WARN]",
    ...,
    "\n"
  )
}


.safe_dir <- function(path) {
  
  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  path
}


.to_logical <- function(
    x,
    default = FALSE
) {
  
  if (
    is.null(x) ||
    length(x) == 0L ||
    is.na(x)
  ) {
    return(default)
  }
  
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  
  tolower(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "true",
    "t",
    "yes",
    "y",
    "1"
  )
}


.is_missing_path <- function(x) {
  
  is.null(x) ||
    length(x) == 0L ||
    is.na(x) ||
    !nzchar(
      trimws(
        as.character(x)
      )
    )
}


.normalize_existing_path <- function(path) {
  
  if (.is_missing_path(path)) {
    return(NULL)
  }
  
  path <- path.expand(
    trimws(
      as.character(path)
    )
  )
  
  normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
}


.safe_scalar_character <- function(
    x,
    default = NA_character_
) {
  
  if (
    is.null(x) ||
    length(x) == 0L ||
    all(is.na(x))
  ) {
    return(default)
  }
  
  value <- trimws(
    as.character(x[[1L]])
  )
  
  if (
    is.na(value) ||
    !nzchar(value)
  ) {
    return(default)
  }
  
  value
}


# ============================================================ #
# Barcode handling
# ============================================================ #

.canonical_barcode <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  # Extract a standard 10X barcode from the end of the string.
  #
  # Examples:
  #   FN_S1256_AAACCCAAGAGCAGCT-1
  #   FN_S1256#AAACCCAAGAGCAGCT-1
  #   FN_S1256:AAACCCAAGAGCAGCT-1
  #
  # become:
  #   AAACCCAAGAGCAGCT-1
  
  extracted <- sub(
    "^.*?([ACGTN]+-[0-9]+)$",
    "\\1",
    x,
    perl = TRUE
  )
  
  is_standard <- grepl(
    "^[ACGTN]+-[0-9]+$",
    extracted
  )
  
  # Fallback for nonstandard barcode formats.
  fallback <- sub(
    "^.*[#_:]",
    "",
    x
  )
  
  output <- ifelse(
    is_standard,
    extracted,
    fallback
  )
  
  trimws(
    as.character(output)
  )
}


.barcode_index <- function(seurat_cells) {
  
  data.frame(
    seurat_cell = seurat_cells,
    barcode = .canonical_barcode(
      seurat_cells
    ),
    stringsAsFactors = FALSE,
    row.names = seurat_cells
  )
}


# ============================================================ #
# Table utilities
# ============================================================ #

.read_table_auto <- function(path) {
  
  path <- .normalize_existing_path(
    path
  )
  
  if (
    is.null(path) ||
    !file.exists(path)
  ) {
    stop(
      "Table not found: ",
      path,
      call. = FALSE
    )
  }
  
  extension <- tolower(
    tools::file_ext(path)
  )
  
  data <- switch(
    extension,
    
    csv = utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    tsv = utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    txt = utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
  
  if (nrow(data) == 0L) {
    stop(
      "Input table is empty: ",
      path,
      call. = FALSE
    )
  }
  
  data
}


.find_column <- function(
    data,
    candidates,
    required = TRUE,
    ignore_case = TRUE
) {
  
  columns <- colnames(data)
  
  if (isTRUE(ignore_case)) {
    
    index <- match(
      tolower(candidates),
      tolower(columns)
    )
    
  } else {
    
    index <- match(
      candidates,
      columns
    )
  }
  
  index <- index[
    !is.na(index)
  ]
  
  if (length(index) > 0L) {
    return(
      columns[
        index[[1L]]
      ]
    )
  }
  
  if (isTRUE(required)) {
    stop(
      "Could not find any of these columns: ",
      paste(
        candidates,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  NULL
}


.clean_assignment <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  missing_assignment <- is.na(x) |
    !nzchar(x) |
    tolower(x) %in% c(
      "na",
      "nan",
      "none",
      "null",
      "unknown",
      "unassigned",
      "negative",
      "neg",
      "not assigned",
      "not_assigned"
    )
  
  x[
    missing_assignment
  ] <- NA_character_
  
  x
}


.is_doublet_label <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  result <- grepl(
    paste(
      c(
        "doublet",
        "multiplet",
        "ambiguous",
        "conflict"
      ),
      collapse = "|"
    ),
    x,
    perl = TRUE
  )
  
  result[
    is.na(result)
  ] <- FALSE
  
  result
}


.is_unassigned_label <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  is.na(x) |
    !nzchar(x) |
    x %in% c(
      "na",
      "nan",
      "none",
      "null",
      "unassigned",
      "unknown",
      "negative",
      "neg",
      "not assigned",
      "not_assigned"
    )
}


# ============================================================ #
# Demuxafy discovery and parsing
# ============================================================ #

.find_demuxafy_table <- function(demuxafy_dir) {
  
  if (.is_missing_path(demuxafy_dir)) {
    return(NULL)
  }
  
  demuxafy_path <- .normalize_existing_path(
    demuxafy_dir
  )
  
  # A direct table path may be supplied instead of a directory.
  if (
    file.exists(demuxafy_path) &&
    !dir.exists(demuxafy_path)
  ) {
    return(demuxafy_path)
  }
  
  if (!dir.exists(demuxafy_path)) {
    return(NULL)
  }
  
  # Explicit locations used by Demuxafy combine_majoritySinglet.
  preferred_paths <- c(
    file.path(
      demuxafy_path,
      "combine_majoritySinglet",
      "combined_results_w_combined_assignments_clean.tsv"
    ),
    file.path(
      demuxafy_path,
      "combine_majoritySinglet",
      "combined_results_w_combined_assignments.tsv"
    ),
    file.path(
      demuxafy_path,
      "combined_results_w_combined_assignments_clean.tsv"
    ),
    file.path(
      demuxafy_path,
      "combined_results_w_combined_assignments.tsv"
    ),
    file.path(
      demuxafy_path,
      "combined_results_with_bulk_mapping.tsv"
    ),
    file.path(
      demuxafy_path,
      "combined_results.tsv"
    )
  )
  
  existing_preferred <- preferred_paths[
    file.exists(
      preferred_paths
    )
  ]
  
  if (length(existing_preferred) > 0L) {
    return(
      existing_preferred[[1L]]
    )
  }
  
  all_tables <- list.files(
    demuxafy_path,
    pattern = "\\.(tsv|txt|csv)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(all_tables) == 0L) {
    return(NULL)
  }
  
  filenames <- tolower(
    basename(all_tables)
  )
  
  # Score candidate tables. Clean combined assignments receive
  # the highest priority.
  score <- integer(
    length(all_tables)
  )
  
  score <- score +
    10L * grepl(
      "combined_results_w_combined_assignments_clean",
      filenames
    ) +
    8L * grepl(
      "combined_results_w_combined_assignments",
      filenames
    ) +
    6L * grepl(
      "combined.*assign",
      filenames
    ) +
    4L * grepl(
      "majority",
      filenames
    ) +
    3L * grepl(
      "assignment",
      filenames
    ) +
    2L * grepl(
      "combined",
      filenames
    ) +
    1L * grepl(
      "result",
      filenames
    )
  
  all_tables[
    order(
      score,
      decreasing = TRUE
    )
  ][[1L]]
}


.detect_barcode_column <- function(data) {
  
  .find_column(
    data,
    candidates = c(
      "Barcode",
      "barcode",
      "BARCODE",
      "cell",
      "Cell",
      "cell_id",
      "cellID",
      "cell_barcode",
      "CellBarcode",
      "Cell_Barcode",
      "CB"
    ),
    required = TRUE
  )
}



.detect_demux_assignment_column <- function(data) {
  
  exact_candidates <- c(
    "MajoritySinglet_Individual_Assignment",
    "majoritysinglet_individual_assignment",
    "demux_status",
    "demux_assignment",
    "combined_assignment",
    "combined.assignments",
    "combined_assignments",
    "final_assignment",
    "final.assignments",
    "majority_assignment",
    "majority.assignments",
    "majoritySinglet",
    "majority_singlet",
    "assignment",
    "best_singlet",
    "singlet_id",
    "donor_id",
    "donor",
    "sample_id",
    "sample"
  )
  
  exact_match <- .find_column(
    data,
    candidates = exact_candidates,
    required = FALSE
  )
  
  if (!is.null(exact_match)) {
    return(exact_match)
  }
  
  columns <- colnames(data)
  
  lower_columns <- tolower(
    columns
  )
  
  index <- grep(
    paste(
      c(
        "majoritysinglet.*individual.*assignment",
        "majority.*individual.*assignment",
        "combined.*assign",
        "majority.*assign",
        "final.*assign",
        "best.*singlet",
        "singlet.*id",
        "donor"
      ),
      collapse = "|"
    ),
    lower_columns
  )
  
  if (length(index) > 0L) {
    return(
      columns[index[[1L]]]
    )
  }
  
  NULL
}


.detect_demux_class_column <- function(data) {
  
  candidates <- c(
    "MajoritySinglet_DropletType",
    "majoritysinglet_droplettype",
    "demux_class",
    "class",
    "classification",
    "combined_classification",
    "combined.classification",
    "doublet_status",
    "doublet.class",
    "doublet_class",
    "call",
    "droplet_type",
    "droplet.class",
    "multiplet_status"
  )
  
  exact_match <- .find_column(
    data,
    candidates = candidates,
    required = FALSE
  )
  
  if (!is.null(exact_match)) {
    return(exact_match)
  }
  
  columns <- colnames(data)
  
  lower_columns <- tolower(
    columns
  )
  
  index <- grep(
    paste(
      c(
        "majoritysinglet.*droplettype",
        "majority.*droplet.*type",
        "combined.*class",
        "doublet.*status",
        "droplet.*type",
        "multiplet"
      ),
      collapse = "|"
    ),
    lower_columns
  )
  
  if (length(index) > 0L) {
    return(
      columns[index[[1L]]]
    )
  }
  
  NULL
}



.prepare_metadata_table <- function(
    data,
    barcode_col,
    prefix = NULL
) {
  
  data <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )
  
  data$barcode <- .canonical_barcode(
    data[[barcode_col]]
  )
  
  data <- data[
    !is.na(data$barcode) &
      nzchar(data$barcode),
    ,
    drop = FALSE
  ]
  
  duplicated_barcodes <- duplicated(
    data$barcode
  )
  
  if (any(duplicated_barcodes)) {
    
    .warn(
      "Removing",
      sum(duplicated_barcodes),
      "duplicated barcode row(s) from assignment table"
    )
    
    data <- data[
      !duplicated_barcodes,
      ,
      drop = FALSE
    ]
  }
  
  if (
    barcode_col %in% colnames(data) &&
    barcode_col != "barcode"
  ) {
    data[[barcode_col]] <- NULL
  }
  
  if (!is.null(prefix)) {
    
    columns_to_prefix <- setdiff(
      colnames(data),
      "barcode"
    )
    
    prefixed_names <- paste0(
      prefix,
      columns_to_prefix
    )
    
    colnames(data)[
      match(
        columns_to_prefix,
        colnames(data)
      )
    ] <- prefixed_names
  }
  
  rownames(data) <- data$barcode
  
  data
}


.add_metadata_by_barcode <- function(
    object,
    metadata,
    overwrite = FALSE
) {
  
  if (!"barcode" %in% colnames(metadata)) {
    stop(
      "Metadata table does not contain a canonical 'barcode' column.",
      call. = FALSE
    )
  }
  
  index <- .barcode_index(
    colnames(object)
  )
  
  matched_index <- match(
    index$barcode,
    metadata$barcode
  )
  
  matched_count <- sum(
    !is.na(matched_index)
  )
  
  if (matched_count == 0L) {
    stop(
      paste0(
        "No barcodes from the assignment table matched ",
        "the Seurat object."
      ),
      call. = FALSE
    )
  }
  
  metadata_columns <- setdiff(
    colnames(metadata),
    "barcode"
  )
  
  aligned_metadata <- metadata[
    matched_index,
    metadata_columns,
    drop = FALSE
  ]
  
  rownames(aligned_metadata) <- index$seurat_cell
  
  for (column in metadata_columns) {
    
    if (
      column %in% colnames(object@meta.data) &&
      !isTRUE(overwrite)
    ) {
      next
    }
    
    object[[column]] <- aligned_metadata[[column]]
  }
  
  list(
    object = object,
    matched = matched_count,
    unmatched = ncol(object) - matched_count
  )
}


# ============================================================ #
# Bulk-assignment column detection
# ============================================================ #

.is_generic_demux_assignment <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  missing <- is.na(x) |
    !nzchar(x)
  
  generic <- grepl(
    "^donor[._ -]*[0-9]+$",
    x,
    perl = TRUE
  ) |
    grepl(
      "^sample[._ -]*[0-9]+$",
      x,
      perl = TRUE
    ) |
    x %in% c(
      "doublet",
      "multiplet",
      "ambiguous",
      "conflict",
      "unassigned",
      "unknown",
      "negative",
      "singlet"
    )
  
  result <- generic & !missing
  
  result[
    is.na(result)
  ] <- FALSE
  
  result
}


.bulk_assignment_quality <- function(x) {
  
  x <- .clean_assignment(
    x
  )
  
  valid <- !is.na(x)
  
  if (!any(valid)) {
    return(
      list(
        score = -Inf,
        valid_fraction = 0,
        generic_fraction = 1,
        biological_fraction = 0
      )
    )
  }
  
  values <- x[valid]
  
  generic <- .is_generic_demux_assignment(
    values
  )
  
  # Biological labels frequently include developmental stage,
  # sample name or donor identifiers such as:
  # S1-PCW7, S2-PCW8, S3-PCW10, FN_PCW6, CTRL-01, etc.
  biological <- grepl(
    paste(
      c(
        "pcw[0-9]",
        "^s[0-9]+[-_]",
        "[-_][a-z]*[0-9]+",
        "ctrl",
        "control",
        "case",
        "patient",
        "sample"
      ),
      collapse = "|"
    ),
    tolower(values),
    perl = TRUE
  ) &
    !generic
  
  valid_fraction <- mean(valid)
  generic_fraction <- mean(generic)
  biological_fraction <- mean(biological)
  
  score <-
    100 * valid_fraction +
    150 * biological_fraction -
    200 * generic_fraction
  
  list(
    score = score,
    valid_fraction = valid_fraction,
    generic_fraction = generic_fraction,
    biological_fraction = biological_fraction
  )
}


.detect_bulk_assignment_column <- function(
    data,
    preferred = NULL,
    verbose = TRUE
) {
  
  log <- .msg(
    verbose
  )
  
  if (
    is.null(data) ||
    !is.data.frame(data) ||
    ncol(data) == 0L
  ) {
    stop(
      "Bulk mapping table is empty or invalid.",
      call. = FALSE
    )
  }
  
  columns <- colnames(
    data
  )
  
  # Explicit column supplied by the user.
  if (
    !is.null(preferred) &&
    length(preferred) > 0L &&
    !is.na(preferred[[1L]]) &&
    nzchar(
      trimws(
        as.character(preferred[[1L]])
      )
    )
  ) {
    
    explicit_column <- .find_column(
      data,
      candidates = as.character(
        preferred[[1L]]
      ),
      required = TRUE,
      ignore_case = TRUE
    )
    
    log(
      "Explicit bulk assignment column:",
      explicit_column
    )
    
    return(
      explicit_column
    )
  }
  
  candidate_names <- c(
    "BulkSample",
    "bulk_sample",
    "bulk.sample",
    "Bulk_Sample",
    "BulkAssignment",
    "bulk_assignment",
    "bulk.assignment",
    "BiologicalSample",
    "biological_sample",
    "biological.sample",
    "MappedSample",
    "mapped_sample",
    "mapped.sample",
    "SampleName",
    "sample_name",
    "sample.name",
    "DonorName",
    "donor_name",
    "donor.name",
    "DonorID",
    "donor_id",
    "donor.id",
    "GenotypeSample",
    "genotype_sample",
    "genotype.sample",
    "Genotype",
    "genotype",
    "sample_id",
    "sample"
  )
  
  matched_columns <- unique(
    columns[
      tolower(columns) %in%
        tolower(candidate_names)
    ]
  )
  
  # Fallback to columns whose names mention biological,
  # bulk, genotype, mapped donor or sample.
  pattern_columns <- columns[
    grepl(
      paste(
        c(
          "bulk",
          "biological.*sample",
          "mapped.*sample",
          "sample.*name",
          "donor.*name",
          "donor.*id",
          "genotype"
        ),
        collapse = "|"
      ),
      tolower(columns),
      perl = TRUE
    )
  ]
  
  matched_columns <- unique(
    c(
      matched_columns,
      pattern_columns
    )
  )
  
  if (length(matched_columns) == 0L) {
    
    stop(
      paste0(
        "Could not identify a bulk biological assignment column. ",
        "Available columns: ",
        paste(
          columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  
  candidate_scores <- vapply(
    matched_columns,
    function(column) {
      
      quality <- .bulk_assignment_quality(
        data[[column]]
      )
      
      name_priority <- match(
        tolower(column),
        tolower(candidate_names)
      )
      
      if (is.na(name_priority)) {
        name_bonus <- 0
      } else {
        name_bonus <- max(
          0,
          50 - name_priority
        )
      }
      
      quality$score + name_bonus
    },
    numeric(1L)
  )
  
  selected_column <- matched_columns[
    which.max(candidate_scores)
  ]
  
  selected_quality <- .bulk_assignment_quality(
    data[[selected_column]]
  )
  
  log(
    "Bulk assignment candidates:",
    paste(
      paste0(
        matched_columns,
        "=",
        round(candidate_scores, 2)
      ),
      collapse = ", "
    )
  )
  
  log(
    "Selected bulk assignment column:",
    selected_column
  )
  
  log(
    "Selected bulk valid fraction:",
    round(
      selected_quality$valid_fraction,
      4
    ),
    "| Generic-label fraction:",
    round(
      selected_quality$generic_fraction,
      4
    ),
    "| Biological-label fraction:",
    round(
      selected_quality$biological_fraction,
      4
    )
  )
  
  selected_column
}


.validate_bulk_assignments <- function(
    x,
    sample,
    column,
    allow_generic = FALSE
) {
  
  assignments <- .clean_assignment(
    x
  )
  
  valid <- !is.na(
    assignments
  )
  
  if (!any(valid)) {
    
    stop(
      paste0(
        "Bulk assignment column '",
        column,
        "' contains no valid assignment for sample ",
        sample,
        "."
      ),
      call. = FALSE
    )
  }
  
  generic <- .is_generic_demux_assignment(
    assignments[valid]
  )
  
  generic_fraction <- mean(
    generic
  )
  
  if (
    !isTRUE(allow_generic) &&
    generic_fraction >= 0.5
  ) {
    
    examples <- unique(
      assignments[valid]
    )
    
    examples <- head(
      examples,
      10L
    )
    
    stop(
      paste0(
        "The selected bulk assignment column '",
        column,
        "' contains mostly generic Demuxafy labels ",
        "rather than biological donor names. Generic fraction: ",
        round(generic_fraction, 3),
        ". Example values: ",
        paste(
          examples,
          collapse = ", "
        ),
        ". Select the column containing labels such as ",
        "S1-PCW7, S2-PCW8, S3-PCW10, etc."
      ),
      call. = FALSE
    )
  }
  
  invisible(
    assignments
  )
}



# ============================================================ #
# Bulk-genotype mapping
# ============================================================ #

.add_bulk_mapping <- function(
    object,
    bulk_tsv,
    sample,
    assignment_column = NULL,
    allow_generic_assignments = FALSE,
    keep_extra_columns = TRUE,
    verbose = TRUE
) {
  
  log <- .msg(
    verbose
  )
  
  bulk <- .read_table_auto(
    bulk_tsv
  )
  
  barcode_col <- .detect_barcode_column(
    bulk
  )
  
  assignment_col <- .detect_bulk_assignment_column(
    data = bulk,
    preferred = assignment_column,
    verbose = verbose
  )
  
  biological_assignments <- .validate_bulk_assignments(
    x = bulk[[assignment_col]],
    sample = sample,
    column = assignment_col,
    allow_generic = allow_generic_assignments
  )
  
  bulk_metadata <- data.frame(
    barcode = .canonical_barcode(
      bulk[[barcode_col]]
    ),
    BulkSample = biological_assignments,
    stringsAsFactors = FALSE
  )
  
  if (isTRUE(keep_extra_columns)) {
    
    extra_columns <- setdiff(
      colnames(bulk),
      c(
        barcode_col,
        assignment_col
      )
    )
    
    if (length(extra_columns) > 0L) {
      
      extra_metadata <- bulk[
        ,
        extra_columns,
        drop = FALSE
      ]
      
      colnames(extra_metadata) <- paste0(
        "bulk_",
        colnames(extra_metadata)
      )
      
      duplicated_names <- duplicated(
        colnames(extra_metadata)
      )
      
      if (any(duplicated_names)) {
        
        extra_metadata <- extra_metadata[
          ,
          !duplicated_names,
          drop = FALSE
        ]
      }
      
      bulk_metadata <- cbind(
        bulk_metadata,
        extra_metadata
      )
    }
  }
  
  valid_barcode <- !is.na(
    bulk_metadata$barcode
  ) &
    nzchar(
      bulk_metadata$barcode
    )
  
  bulk_metadata <- bulk_metadata[
    valid_barcode,
    ,
    drop = FALSE
  ]
  
  duplicated_barcode <- duplicated(
    bulk_metadata$barcode
  )
  
  if (any(duplicated_barcode)) {
    
    .warn(
      "Removing",
      sum(duplicated_barcode),
      "duplicated barcode row(s) from bulk mapping"
    )
    
    bulk_metadata <- bulk_metadata[
      !duplicated_barcode,
      ,
      drop = FALSE
    ]
  }
  
  added <- .add_metadata_by_barcode(
    object = object,
    metadata = bulk_metadata,
    overwrite = TRUE
  )
  
  log(
    "Bulk assignment column:",
    assignment_col
  )
  
  log(
    "Bulk mapping matched:",
    added$matched,
    "| Unmatched:",
    added$unmatched
  )
  
  matched_assignments <- added$object$BulkSample
  
  log(
    "Bulk assignment examples:",
    paste(
      head(
        unique(
          stats::na.omit(
            matched_assignments
          )
        ),
        10L
      ),
      collapse = ", "
    )
  )
  
  list(
    object = added$object,
    matched = added$matched,
    unmatched = added$unmatched,
    assignment_column = assignment_col
  )
}


# ============================================================ #
# Main reusable function
# ============================================================ #

run_demux <- function(
    sample,
    rds_in,
    demuxafy_dir = NULL,
    bulk_tsv = NULL,
    remove_doublets = FALSE,
    remove_unassigned = TRUE,
    assignment_preference = c(
      "demux",
      "demuxafy",
      "bulk"
    ),
    bulk_assignment_column = NULL,
    allow_generic_bulk_assignments = FALSE,
    keep_original_metadata = TRUE,
    metadata_mode = c(
      "full",
      "minimal"
    ),
    outdir = "results/02_demuxed",
    verbose = TRUE
) {
  
  log <- .msg(
    verbose
  )
  
  # -------------------------------------------------------- #-- #
  # Validate and normalize parameters
  # -------------------------------------------------------- #-- #
  
  sample <- .safe_scalar_character(
    sample
  )
  
  if (is.na(sample)) {
    
    stop(
      "A non-empty sample name is required.",
      call. = FALSE
    )
  }
  
  assignment_preference <- match.arg(
    assignment_preference
  )
  
  metadata_mode <- match.arg(
    metadata_mode
  )
  
  if (!isTRUE(keep_original_metadata)) {
    metadata_mode <- "minimal"
  }
  
  if (identical(
    assignment_preference,
    "demuxafy"
  )) {
    assignment_preference <- "demux"
  }
  
  remove_doublets <- .to_logical(
    remove_doublets,
    default = FALSE
  )
  
  remove_unassigned <- .to_logical(
    remove_unassigned,
    default = TRUE
  )
  
  allow_generic_bulk_assignments <- .to_logical(
    allow_generic_bulk_assignments,
    default = FALSE
  )
  
  rds_in <- .normalize_existing_path(
    rds_in
  )
  
  demuxafy_dir <- .normalize_existing_path(
    demuxafy_dir
  )
  
  bulk_tsv <- .normalize_existing_path(
    bulk_tsv
  )
  
  outdir <- .safe_dir(
    outdir
  )
  
  if (
    is.null(rds_in) ||
    !file.exists(rds_in)
  ) {
    
    stop(
      "Input RDS not found: ",
      rds_in,
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------- #-- #
  # Read input object
  # -------------------------------------------------------- #-- #
  
  log(
    "==================================="
  )
  
  log(
    "DEMULTIPLEXING METADATA"
  )
  
  log(
    "==================================="
  )
  
  log(
    "Sample:",
    sample
  )
  
  log(
    "Input:",
    rds_in
  )
  
  log(
    "Assignment preference:",
    assignment_preference
  )
  
  object <- readRDS(
    rds_in
  )
  
  if (!inherits(
    object,
    "Seurat"
  )) {
    
    stop(
      "Input RDS does not contain a Seurat object.",
      call. = FALSE
    )
  }
  
  if (ncol(object) == 0L) {
    
    stop(
      "Input Seurat object contains no cells.",
      call. = FALSE
    )
  }
  
  object$sample <- sample
  object$orig.ident <- sample
  
  initial_cells <- ncol(
    object
  )
  
  object$original_cell_name <- colnames(
    object
  )
  
  object$barcode <- .canonical_barcode(
    colnames(object)
  )
  
  # -------------------------------------------------------- #-- #
  # Discover assignment inputs
  # -------------------------------------------------------- #-- #
  
  demux_table <- .find_demuxafy_table(
    demuxafy_dir
  )
  
  has_demux_table <- !is.null(
    demux_table
  ) &&
    file.exists(
      demux_table
    )
  
  has_bulk_table <- !is.null(
    bulk_tsv
  ) &&
    file.exists(
      bulk_tsv
    )
  
  if (
    !has_demux_table &&
    !has_bulk_table
  ) {
    
    stop(
      paste0(
        "No valid Demuxafy assignment table or bulk-mapping ",
        "table was found for sample: ",
        sample,
        ". Demuxafy directory: ",
        demuxafy_dir %||%
          "<not supplied>",
        "; bulk table: ",
        bulk_tsv %||%
          "<not supplied>"
      ),
      call. = FALSE
    )
  }
  
  if (
    identical(
      assignment_preference,
      "bulk"
    ) &&
    !has_bulk_table
  ) {
    
    stop(
      paste0(
        "assignment_preference='bulk', but no valid bulk ",
        "mapping table was found for sample: ",
        sample
      ),
      call. = FALSE
    )
  }
  
  demux_added <- FALSE
  bulk_added <- FALSE
  
  demux_matched <- 0L
  demux_unmatched <- initial_cells
  
  bulk_matched <- 0L
  bulk_unmatched <- initial_cells
  selected_bulk_column <- NA_character_
  
  
  # -------------------------------------------------------- #-- #
  # Add Demuxafy metadata
  # -------------------------------------------------------- #-- #
  
  if (has_demux_table) {
    
    demux_table <- normalizePath(
      demux_table,
      winslash = "/",
      mustWork = TRUE
    )
    
    log(
      "Demuxafy table:",
      demux_table
    )
    
    demux <- .read_table_auto(
      demux_table
    )
    
    barcode_col <- .detect_barcode_column(
      demux
    )
    
    assignment_col <- .detect_demux_assignment_column(
      demux
    )
    
    class_col <- .detect_demux_class_column(
      demux
    )
    
    demux_metadata <- .prepare_metadata_table(
      data = demux,
      barcode_col = barcode_col,
      prefix = "demux_"
    )
    
    # -------------------------------------------------------- #
    # Standardized Demuxafy assignment
    # -------------------------------------------------------- #
    
    if (!is.null(assignment_col)) {
      
      prefixed_assignment_col <- paste0(
        "demux_",
        assignment_col
      )
      
      if (
        prefixed_assignment_col %in%
        colnames(demux_metadata)
      ) {
        
        demux_metadata$demux_status <-
          .clean_assignment(
            demux_metadata[[prefixed_assignment_col]]
          )
        
      } else {
        
        demux_metadata$demux_status <-
          NA_character_
        
        .warn(
          "Detected assignment column was not retained:",
          assignment_col
        )
      }
      
    } else {
      
      demux_metadata$demux_status <-
        NA_character_
      
      .warn(
        paste0(
          "No clear Demuxafy assignment column was detected. ",
          "Original Demuxafy columns were retained."
        )
      )
    }
    
    # -------------------------------------------------------- #
    # Standardized Demuxafy classification
    # -------------------------------------------------------- #
    
    if (!is.null(class_col)) {
      
      prefixed_class_col <- paste0(
        "demux_",
        class_col
      )
      
      if (
        prefixed_class_col %in%
        colnames(demux_metadata)
      ) {
        
        demux_metadata$demux_class <-
          as.character(
            demux_metadata[[prefixed_class_col]]
          )
        
      } else {
        
        demux_metadata$demux_class <-
          NA_character_
        
        .warn(
          "Detected Demuxafy class column was not retained:",
          class_col
        )
      }
      
    } else {
      
      demux_metadata$demux_class <-
        NA_character_
    }
    
    # -------------------------------------------------------- #
    # Add aligned Demuxafy metadata to Seurat object
    # -------------------------------------------------------- #
    
    added <- .add_metadata_by_barcode(
      object = object,
      metadata = demux_metadata,
      overwrite = FALSE
    )
    
    object <- added$object
    
    demux_matched <- added$matched
    demux_unmatched <- added$unmatched
    demux_added <- TRUE
    
    log(
      "Demuxafy assignment column:",
      assignment_col %||%
        "<not detected>"
    )
    
    log(
      "Demuxafy class column:",
      class_col %||%
        "<not detected>"
    )
    
    log(
      "Demuxafy matched:",
      demux_matched,
      "| Unmatched:",
      demux_unmatched
    )
    
  } else if (!is.null(demuxafy_dir)) {
    
    .warn(
      "No Demuxafy result table found under:",
      demuxafy_dir
    )
  }
  
  
  # -------------------------------------------------------- #-- #
  # Add bulk-genotype metadata
  # -------------------------------------------------------- #-- #
  
  if (has_bulk_table) {
    
    bulk_tsv <- normalizePath(
      bulk_tsv,
      winslash = "/",
      mustWork = TRUE
    )
    
    log(
      "Bulk mapping:",
      bulk_tsv
    )
    
    bulk_result <- .add_bulk_mapping(
      object = object,
      bulk_tsv = bulk_tsv,
      sample = sample,
      assignment_column = bulk_assignment_column,
      allow_generic_assignments =
        allow_generic_bulk_assignments,
      keep_extra_columns =
        identical(
          metadata_mode,
          "full"
        ),
      verbose = verbose
    )
    
    object <- bulk_result$object
    
    bulk_matched <- bulk_result$matched
    bulk_unmatched <- bulk_result$unmatched
    selected_bulk_column <-
      bulk_result$assignment_column
    
    bulk_added <- TRUE
    
  } else if (!is.null(bulk_tsv)) {
    
    .warn(
      "Bulk mapping table not found:",
      bulk_tsv
    )
  }
  
  # -------------------------------------------------------- #-- #
  # Construct standardized final assignment
  # -------------------------------------------------------- #-- #
  
  bulk_status <- if (
    "BulkSample" %in% colnames(object@meta.data)
  ) {
    
    .clean_assignment(
      object@meta.data[["BulkSample"]]
    )
    
  } else {
    
    rep(
      NA_character_,
      ncol(object)
    )
  }
  
  demux_status <- if (
    "demux_status" %in% colnames(object@meta.data)
  ) {
    
    .clean_assignment(
      object@meta.data[["demux_status"]]
    )
    
  } else {
    
    rep(
      NA_character_,
      ncol(object)
    )
  }
  
  log(
    "Requested assignment preference:",
    assignment_preference
  )
  
  log(
    "BulkSample available:",
    "BulkSample" %in% colnames(object@meta.data)
  )
  
  log(
    "Bulk assignment examples before final selection:",
    paste(
      head(
        unique(
          stats::na.omit(
            bulk_status
          )
        ),
        10L
      ),
      collapse = ", "
    )
  )
  
  log(
    "Demux assignment examples before final selection:",
    paste(
      head(
        unique(
          stats::na.omit(
            demux_status
          )
        ),
        10L
      ),
      collapse = ", "
    )
  )
  
  if (identical(
    assignment_preference,
    "bulk"
  )) {
    
    if (!bulk_added) {
      
      stop(
        paste0(
          "assignment_preference='bulk', but bulk metadata ",
          "was not added for sample ",
          sample,
          "."
        ),
        call. = FALSE
      )
    }
    
    final_status <- bulk_status
    
  } else {
    
    final_status <- demux_status
    
    # Demux mode may use the biological bulk label only when
    # the Demuxafy assignment is missing.
    fill_from_bulk <- is.na(
      final_status
    ) &
      !is.na(
        bulk_status
      )
    
    final_status[
      fill_from_bulk
    ] <- bulk_status[
      fill_from_bulk
    ]
  }
  
  final_status <- .clean_assignment(
    final_status
  )
  
  # In bulk mode, never replace missing biological assignments
  # with donor0, donor1, donor2, etc.
  if (identical(
    assignment_preference,
    "bulk"
  )) {
    
    generic_bulk <- !is.na(final_status) &
      .is_generic_demux_assignment(
        final_status
      )
    
    if (any(generic_bulk)) {
      
      bad_examples <- unique(
        final_status[
          generic_bulk
        ]
      )
      
      stop(
        paste0(
          "Bulk mode produced generic donor labels: ",
          paste(
            head(
              bad_examples,
              10L
            ),
            collapse = ", "
          ),
          ". Expected biological labels from BulkSample, ",
          "such as S1-PCW7 or S6-PCW12."
        ),
        call. = FALSE
      )
    }
  }
  
  object$status <- final_status
  object$donor_id <- final_status
  
  log(
    "Final assignment source:",
    if (
      identical(
        assignment_preference,
        "bulk"
      )
    ) {
      "BulkSample"
    } else {
      "demux_status"
    }
  )
  
  log(
    "Final assignment examples:",
    paste(
      head(
        unique(
          stats::na.omit(
            object$status
          )
        ),
        10L
      ),
      collapse = ", "
    )
  )
  
  # -------------------------------------------------------- #-- #
  # Identify doublets and unassigned cells
  # -------------------------------------------------------- #-- #
  
  is_doublet <- rep(
    FALSE,
    ncol(object)
  )
  
  candidate_doublet_columns <- intersect(
    c(
      "demux_class",
      "demux_status",
      "BulkSample",
      "status",
      "demux_combined_classification",
      "demux_combined.classification",
      "demux_doublet_status",
      "demux_doublet.class",
      "demux_classification"
    ),
    colnames(object@meta.data)
  )
  
  for (column in candidate_doublet_columns) {
    
    is_doublet <- is_doublet |
      .is_doublet_label(
        object@meta.data[[column]]
      )
  }
  
  object$demux_is_doublet <- is_doublet
  
  object$demux_is_unassigned <-
    .is_unassigned_label(
      object$status
    )
  
  # -------------------------------------------------------- #-- #
  # Filter cells
  # -------------------------------------------------------- #-- #
  
  keep <- rep(
    TRUE,
    ncol(object)
  )
  
  if (isTRUE(remove_doublets)) {
    
    keep <- keep &
      !object$demux_is_doublet
  }
  
  if (isTRUE(remove_unassigned)) {
    
    keep <- keep &
      !object$demux_is_unassigned
  }
  
  keep[
    is.na(keep)
  ] <- FALSE
  
  removed_doublets <- sum(
    object$demux_is_doublet &
      !keep,
    na.rm = TRUE
  )
  
  removed_unassigned <- sum(
    object$demux_is_unassigned &
      !keep,
    na.rm = TRUE
  )
  
  kept_is_doublet <- as.logical(
    object$demux_is_doublet[
      keep
    ]
  )
  
  kept_is_unassigned <- as.logical(
    object$demux_is_unassigned[
      keep
    ]
  )
  
  object <- object[
    ,
    keep,
    drop = FALSE
  ]
  
  final_cells <- ncol(
    object
  )
  
  if (final_cells == 0L) {
    
    stop(
      paste0(
        "No cells remained after demultiplexing filters for ",
        "sample: ",
        sample
      ),
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------- #-- #
  # Standardize metadata
  # -------------------------------------------------------- #-- #
  
  if (
    !"species" %in%
    colnames(object@meta.data)
  ) {
    
    object$species <- rep(
      NA_character_,
      ncol(object)
    )
  }
  
  object$status <- .clean_assignment(
    object$status
  )
  
  object$donor_id <- object$status
  
  # -------------------------------------------------------- #-- #
  # Calculate or harmonize mitochondrial percentage
  # -------------------------------------------------------- #-- #
  
  if (
    !"percent_mito" %in%
    colnames(object@meta.data)
  ) {
    
    mitochondrial_candidates <- c(
      "percent.mt",
      "percent.mito",
      "pct_counts_mt"
    )
    
    mitochondrial_column <- intersect(
      mitochondrial_candidates,
      colnames(object@meta.data)
    )
    
    if (length(
      mitochondrial_column
    ) > 0L) {
      
      mitochondrial_column <-
        mitochondrial_column[[1L]]
      
      mitochondrial_values <- as.numeric(
        object@meta.data[[mitochondrial_column]]
      )
      
      # Seurat PercentageFeatureSet usually produces 0–100.
      # Convert to proportions to match FN_PCW6 metadata.
      if (
        any(
          mitochondrial_values > 1,
          na.rm = TRUE
        )
      ) {
        
        mitochondrial_values <-
          mitochondrial_values / 100
      }
      
      object$percent_mito <-
        mitochondrial_values
      
    } else {
      
      assay_name <- if (
        "RNA" %in% names(
          object@assays
        )
      ) {
        "RNA"
      } else {
        SeuratObject::DefaultAssay(
          object
        )
      }
      
      feature_names <- rownames(
        object[[assay_name]]
      )
      
      mitochondrial_features <- grep(
        pattern = "^MT-",
        x = feature_names,
        value = TRUE,
        ignore.case = TRUE
      )
      
      if (
        length(
          mitochondrial_features
        ) > 0L
      ) {
        
        object[["percent_mito"]] <-
          Seurat::PercentageFeatureSet(
            object = object,
            features = mitochondrial_features,
            assay = assay_name
          ) / 100
        
      } else {
        
        .warn(
          "No mitochondrial features matching '^MT-' were found."
        )
        
        object$percent_mito <- rep(
          NA_real_,
          ncol(object)
        )
      }
    }
  }
  
  # -------------------------------------------------------- #-- #
  # Reduce metadata when requested
  # -------------------------------------------------------- #-- #
  
  if (identical(
    metadata_mode,
    "minimal"
  )) {
    
    minimal_columns <- c(
      "orig.ident",
      "nCount_RNA",
      "nFeature_RNA",
      "sample",
      "species",
      "status",
      "donor_id",
      "percent_mito"
    )
    
    missing_columns <- setdiff(
      minimal_columns,
      colnames(object@meta.data)
    )
    
    for (column in missing_columns) {
      
      if (identical(
        column,
        "percent_mito"
      )) {
        
        object@meta.data[[column]] <- rep(
          NA_real_,
          ncol(object)
        )
        
      } else {
        
        object@meta.data[[column]] <- rep(
          NA_character_,
          ncol(object)
        )
      }
    }
    
    object@meta.data <- object@meta.data[
      ,
      minimal_columns,
      drop = FALSE
    ]
  }
  
  log(
    "Metadata mode:",
    metadata_mode
  )
  
  log(
    "Metadata columns:",
    paste(
      colnames(object@meta.data),
      collapse = ", "
    )
  )
  
  # -------------------------------------------------------- #-- #
  # Save output object
  # -------------------------------------------------------- #-- #
  
  rds_out <- file.path(
    outdir,
    paste0(
      sample,
      "_demuxed.rds"
    )
  )
  
  saveRDS(
    object,
    rds_out
  )
  
  # -------------------------------------------------------- #-- #
  # Save run summary
  # -------------------------------------------------------- #-- #
  
  summary <- data.frame(
    Sample = sample,
    Initial_Cells = initial_cells,
    Final_Cells = final_cells,
    Removed_Cells =
      initial_cells - final_cells,
    Removed_Doublets =
      removed_doublets,
    Removed_Unassigned =
      removed_unassigned,
    Demuxafy_Added =
      demux_added,
    Demuxafy_Matched =
      demux_matched,
    Demuxafy_Unmatched =
      demux_unmatched,
    Bulk_Mapping_Added =
      bulk_added,
    Bulk_Mapping_Matched =
      bulk_matched,
    Bulk_Mapping_Unmatched =
      bulk_unmatched,
    Bulk_Assignment_Column =
      selected_bulk_column,
    Assignment_Preference =
      assignment_preference,
    Metadata_Mode =
      metadata_mode,
    Keep_Original_Metadata =
      keep_original_metadata,
    Remove_Doublets =
      remove_doublets,
    Remove_Unassigned =
      remove_unassigned,
    Demuxafy_Table = if (
      has_demux_table
    ) {
      demux_table
    } else {
      NA_character_
    },
    Bulk_Table = if (
      has_bulk_table
    ) {
      bulk_tsv
    } else {
      NA_character_
    },
    Output = rds_out,
    stringsAsFactors = FALSE
  )
  
  summary_file <- file.path(
    outdir,
    paste0(
      sample,
      "_demux_summary.tsv"
    )
  )
  
  utils::write.table(
    summary,
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  
  # -------------------------------------------------------- #-- #
  # Save assignment counts
  # -------------------------------------------------------- #-- #
  
  assignment_summary <- as.data.frame(
    table(
      status = object$status,
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(
    assignment_summary
  ) <- c(
    "status",
    "cells"
  )
  
  assignment_counts_file <- file.path(
    outdir,
    paste0(
      sample,
      "_assignment_counts.tsv"
    )
  )
  
  utils::write.table(
    assignment_summary,
    assignment_counts_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  
  # -------------------------------------------------------- #-- #
  # Save classification counts
  # -------------------------------------------------------- #-- #
  
  classification_summary <- data.frame(
    category = c(
      "assigned_singlet",
      "doublet",
      "unassigned"
    ),
    cells = c(
      sum(
        !kept_is_doublet &
          !kept_is_unassigned,
        na.rm = TRUE
      ),
      sum(
        kept_is_doublet,
        na.rm = TRUE
      ),
      sum(
        kept_is_unassigned,
        na.rm = TRUE
      )
    ),
    stringsAsFactors = FALSE
  )
  
  utils::write.table(
    classification_summary,
    file.path(
      outdir,
      paste0(
        sample,
        "_demux_classification_counts.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # -------------------------------------------------------- #-- #
  # Final messages
  # -------------------------------------------------------- #-- #
  
  log(
    "Cells before demultiplexing:",
    initial_cells
  )
  
  log(
    "Cells after demultiplexing:",
    final_cells
  )
  
  log(
    "Removed doublets:",
    removed_doublets
  )
  
  log(
    "Removed unassigned:",
    removed_unassigned
  )
  
  log(
    "Assignment preference:",
    assignment_preference
  )
  
  log(
    "Bulk assignment column:",
    selected_bulk_column
  )
  
  log(
    "Saved:",
    rds_out
  )
  
  invisible(
    object
  )
}

# ============================================================ #
# Optional Snakemake entry point
# ============================================================ #

# This block is not triggered by:
#
#   source("scripts/02_add_demux.R")
#
# unless a Snakemake object was actually injected into the
# script environment.

if (
  exists(
    "snakemake",
    inherits = FALSE
  )
) {
  
  run_demux(
    sample = snakemake@params[["sample"]],
    rds_in = snakemake@input[["rds"]],
    demuxafy_dir = snakemake@params[["demuxafy_dir"]] %||% NULL,
    bulk_tsv = snakemake@params[["bulk_tsv"]] %||% NULL,
    remove_doublets = snakemake@params[["remove_doublets"]] %||% FALSE,
    remove_unassigned = snakemake@params[["remove_unassigned"]] %||% TRUE,
    assignment_preference = snakemake@params[["assignment_preference"]] %||% "bulk",
    keep_original_metadata = snakemake@params[["keep_original_metadata"]] %||% FALSE,
    metadata_mode = snakemake@params[["metadata_mode"]] %||% "minimal",
    outdir = dirname(snakemake@output[["rds"]])
  )
}