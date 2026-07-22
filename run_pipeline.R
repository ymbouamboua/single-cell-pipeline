#!/usr/bin/env Rscript

# ============================================================ #
# 10x single-cell pipeline standalone runner
#
# Examples:
#
# Rscript run_pipeline.R --workflow basic
# Rscript run_pipeline.R --workflow demux
# Rscript run_pipeline.R --workflow contam
# Rscript run_pipeline.R --workflow full
#
# Rscript run_pipeline.R \
#   --workflow demux \
#   --sample FN_S1256
#
# Rscript run_pipeline.R \
#   --step qc \
#   --sample FN_S1256
#
# Rscript run_pipeline.R \
#   --step integrate \
#   --input_rds results/03_qc/S1_qc.rds,results/03_qc/S2_qc.rds
# ============================================================ #


suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})


# ============================================================ #
# General helpers
# ============================================================ #

`%||%` <- function(x, y) {
  
  if (
    is.null(x) ||
    length(x) == 0L ||
    all(is.na(x))
  ) {
    return(y)
  }
  
  x
}


msg <- function(...) {
  
  cat(
    "[INFO]",
    ...,
    "\n"
  )
}


warn <- function(...) {
  
  cat(
    "[WARN]",
    ...,
    "\n"
  )
}


fail <- function(...) {
  
  stop(
    paste(...),
    call. = FALSE
  )
}


to_bool <- function(
    x,
    default = FALSE
) {
  
  if (
    is.null(x) ||
    length(x) == 0L ||
    all(is.na(x))
  ) {
    return(default)
  }
  
  if (is.logical(x)) {
    return(
      isTRUE(x[[1L]])
    )
  }
  
  value <- tolower(
    trimws(
      as.character(x[[1L]])
    )
  )
  
  if (
    value %in% c(
      "true",
      "t",
      "yes",
      "y",
      "1"
    )
  ) {
    return(TRUE)
  }
  
  if (
    value %in% c(
      "false",
      "f",
      "no",
      "n",
      "0",
      ""
    )
  ) {
    return(FALSE)
  }
  
  default
}


get_col <- function(
    data,
    column,
    default = NA
) {
  
  if (
    !column %in% colnames(data) ||
    nrow(data) == 0L
  ) {
    return(default)
  }
  
  value <- data[[column]][[1L]]
  
  if (
    is.null(value) ||
    length(value) == 0L ||
    is.na(value)
  ) {
    return(default)
  }
  
  if (
    is.character(value) &&
    !nzchar(trimws(value))
  ) {
    return(default)
  }
  
  value
}


is_missing_value <- function(x) {
  
  is.null(x) ||
    length(x) == 0L ||
    all(is.na(x)) ||
    (
      is.character(x) &&
        !nzchar(
          trimws(
            as.character(x[[1L]])
          )
        )
    )
}


build_steps <- function(workflow) {
  
  switch(
    tolower(workflow),
    
    basic = c(
      "load",
      "qc",
      "integrate"
    ),
    
    demux = c(
      "load",
      "demux",
      "qc",
      "integrate"
    ),
    
    contam = c(
      "load",
      "demux",
      "qc",
      "contam",
      "integrate"
    ),
    
    full = c(
      "load",
      "demux",
      "qc",
      "contam",
      "integrate",
      "deg"
    ),
    
    fail(
      "Unknown workflow:",
      workflow
    )
  )
}


split_paths <- function(x) {
  
  paths <- trimws(
    strsplit(
      x,
      ",",
      fixed = TRUE
    )[[1L]]
  )
  
  paths[
    nzchar(paths)
  ]
}


# ============================================================ #
# Detect pipeline repository
# ============================================================ #

get_script_path <- function() {
  
  arguments <- commandArgs(
    trailingOnly = FALSE
  )
  
  file_argument <- grep(
    "^--file=",
    arguments,
    value = TRUE
  )
  
  if (length(file_argument) == 0L) {
    
    return(
      normalizePath(
        ".",
        mustWork = TRUE
      )
    )
  }
  
  script_file <- sub(
    "^--file=",
    "",
    file_argument[[1L]]
  )
  
  dirname(
    normalizePath(
      script_file,
      mustWork = TRUE
    )
  )
}


PIPELINE_DIR <- get_script_path()

candidate_script_dirs <- c(
  file.path(
    PIPELINE_DIR,
    "workflow",
    "scripts"
  ),
  file.path(
    PIPELINE_DIR,
    "scripts"
  )
)

existing_script_dirs <- candidate_script_dirs[
  dir.exists(candidate_script_dirs)
]

if (length(existing_script_dirs) == 0L) {
  
  fail(
    "No pipeline scripts directory found. Checked:",
    paste(
      candidate_script_dirs,
      collapse = ", "
    )
  )
}

SCRIPT_DIR <- existing_script_dirs[[1L]]


source_pipeline_script <- function(filename) {
  
  script_path <- file.path(
    SCRIPT_DIR,
    filename
  )
  
  if (!file.exists(script_path)) {
    
    fail(
      "Pipeline script not found:",
      script_path
    )
  }
  
  source(
    script_path,
    local = globalenv()
  )
}


# ============================================================ #
# Command-line options
# ============================================================ #

option_list <- list(
  
  make_option(
    "--config",
    type = "character",
    default = "config/pipeline_config.yaml",
    help = paste(
      "Project configuration file",
      "[default: %default]"
    )
  ),
  
  make_option(
    "--workflow",
    type = "character",
    default = NULL,
    help = paste(
      "Workflow:",
      "basic, demux, contam or full"
    )
  ),
  
  make_option(
    "--step",
    type = "character",
    default = NULL,
    help = paste(
      "Single step:",
      "load, demux, qc, contam, integrate or deg"
    )
  ),
  
  make_option(
    "--sample",
    type = "character",
    default = NULL,
    help = "Comma-separated SampleID values"
  ),
  
  make_option(
    "--input_rds",
    type = "character",
    default = NULL,
    help = paste(
      "Comma-separated RDS files for integration,",
      "or one integrated RDS for DEG"
    )
  ),
  
  make_option(
    "--continue_on_error",
    type = "logical",
    default = FALSE,
    help = paste(
      "Continue with other samples after",
      "a per-sample failure [default: %default]"
    )
  )
)


options <- parse_args(
  OptionParser(
    option_list = option_list
  )
)


# ============================================================ #
# Configuration
# ============================================================ #

config_file <- normalizePath(
  options$config,
  mustWork = FALSE
)

if (!file.exists(config_file)) {
  
  fail(
    "Configuration file not found:",
    config_file
  )
}

cfg <- yaml::read_yaml(
  config_file
)

if (is.null(cfg$master_summary)) {
  
  fail(
    "Missing required configuration field:",
    "master_summary"
  )
}


# ============================================================ #
# Project-relative path handling
# ============================================================ #

PROJECT_DIR <- dirname(
  dirname(config_file)
)


resolve_project_path <- function(path) {
  
  if (is_missing_value(path)) {
    return(NULL)
  }
  
  path <- path.expand(
    trimws(
      as.character(path[[1L]])
    )
  )
  
  is_absolute <- grepl(
    "^/",
    path
  ) ||
    grepl(
      "^[A-Za-z]:[/\\\\]",
      path
    )
  
  if (!is_absolute) {
    
    path <- file.path(
      PROJECT_DIR,
      path
    )
  }
  
  normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
}


resolve_required_path <- function(
    path,
    description
) {
  
  resolved <- resolve_project_path(
    path
  )
  
  if (
    is.null(resolved) ||
    !file.exists(resolved)
  ) {
    
    fail(
      paste0(
        description,
        " not found:"
      ),
      resolved %||%
        "<not supplied>"
    )
  }
  
  resolved
}


resolve_raw_matrix <- function(cellranger_dir) {
  
  raw_directory <- file.path(
    cellranger_dir,
    "raw_feature_bc_matrix"
  )
  
  raw_h5 <- file.path(
    cellranger_dir,
    "raw_feature_bc_matrix.h5"
  )
  
  if (dir.exists(raw_directory)) {
    return(raw_directory)
  }
  
  if (file.exists(raw_h5)) {
    return(raw_h5)
  }
  
  fail(
    "Raw 10x matrix not found under:",
    cellranger_dir
  )
}


# ============================================================ #
# Master summary
# ============================================================ #

master_file <- resolve_required_path(
  cfg$master_summary,
  "Master summary"
)

master <- read.csv(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "null",
    "NULL"
  )
)


required_master_columns <- c(
  "SampleID",
  "CellRanger_Dir"
)

missing_master_columns <- setdiff(
  required_master_columns,
  colnames(master)
)

if (length(missing_master_columns) > 0L) {
  
  fail(
    "Missing master-summary column(s):",
    paste(
      missing_master_columns,
      collapse = ", "
    )
  )
}

master$SampleID <- trimws(
  as.character(
    master$SampleID
  )
)

if (
  anyNA(master$SampleID) ||
  any(!nzchar(master$SampleID))
) {
  
  fail(
    "Missing or empty SampleID values detected."
  )
}

duplicated_samples <- unique(
  master$SampleID[
    duplicated(master$SampleID)
  ]
)

if (length(duplicated_samples) > 0L) {
  
  fail(
    "Duplicated SampleID values:",
    paste(
      duplicated_samples,
      collapse = ", "
    )
  )
}

if ("Include" %in% colnames(master)) {
  
  include_sample <- vapply(
    master$Include,
    to_bool,
    logical(1L),
    default = TRUE
  )
  
  master <- master[
    include_sample,
    ,
    drop = FALSE
  ]
}

if (nrow(master) == 0L) {
  
  fail(
    "No enabled samples remain in the master summary."
  )
}


# ============================================================ #
# Workflow setup
# ============================================================ #

outdir <- resolve_project_path(
  cfg$outdir %||%
    "results"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

workflow <- options$workflow %||%
  cfg$workflow %||%
  "basic"

steps <- if (!is.null(options$step)) {
  
  tolower(
    trimws(
      options$step
    )
  )
  
} else {
  
  build_steps(
    workflow
  )
}

valid_steps <- c(
  "load",
  "demux",
  "qc",
  "contam",
  "integrate",
  "deg"
)

invalid_steps <- setdiff(
  steps,
  valid_steps
)

if (length(invalid_steps) > 0L) {
  
  fail(
    "Invalid step(s):",
    paste(
      invalid_steps,
      collapse = ", "
    )
  )
}

samples <- if (!is.null(options$sample)) {
  
  split_paths(
    options$sample
  )
  
} else {
  
  master$SampleID
}

unknown_samples <- setdiff(
  samples,
  master$SampleID
)

if (length(unknown_samples) > 0L) {
  
  fail(
    "Unknown sample(s):",
    paste(
      unknown_samples,
      collapse = ", "
    )
  )
}


msg(
  "Pipeline directory:",
  PIPELINE_DIR
)

msg(
  "Scripts directory: ",
  SCRIPT_DIR
)

msg(
  "Project directory: ",
  PROJECT_DIR
)

msg(
  "Configuration:     ",
  config_file
)

msg(
  "Master summary:    ",
  master_file
)

msg(
  "Workflow:          ",
  workflow
)

msg(
  "Steps:             ",
  paste(
    steps,
    collapse = " -> "
  )
)

msg(
  "Samples:           ",
  paste(
    samples,
    collapse = ", "
  )
)


utils_file <- file.path(
  SCRIPT_DIR,
  "utils.R"
)

if (file.exists(utils_file)) {
  
  source(
    utils_file,
    local = globalenv()
  )
}


# ============================================================ #
# Demultiplexing helpers
# ============================================================ #

sample_has_demux <- function(meta) {
  
  demuxafy_value <- get_col(
    meta,
    "Demuxafy_Dir",
    NA_character_
  )
  
  bulk_value <- get_col(
    meta,
    "BulkMapping_TSV",
    NA_character_
  )
  
  !is_missing_value(demuxafy_value) ||
    !is_missing_value(bulk_value)
}


resolve_demux_inputs <- function(
    sample,
    meta
) {
  
  raw_demuxafy <- get_col(
    meta,
    "Demuxafy_Dir",
    NA_character_
  )
  
  raw_bulk <- get_col(
    meta,
    "BulkMapping_TSV",
    NA_character_
  )
  
  demuxafy_dir <- resolve_project_path(
    raw_demuxafy
  )
  
  bulk_tsv <- resolve_project_path(
    raw_bulk
  )
  
  if (
    !is.null(demuxafy_dir) &&
    !dir.exists(demuxafy_dir) &&
    !file.exists(demuxafy_dir)
  ) {
    
    warn(
      sample,
      ":: Demuxafy input not found; ignoring:",
      demuxafy_dir
    )
    
    demuxafy_dir <- NULL
  }
  
  if (
    !is.null(bulk_tsv) &&
    !file.exists(bulk_tsv)
  ) {
    
    warn(
      sample,
      ":: bulk-mapping table not found;",
      "continuing without it:",
      bulk_tsv
    )
    
    bulk_tsv <- NULL
  }
  
  list(
    demuxafy_dir = demuxafy_dir,
    bulk_tsv = bulk_tsv
  )
}


get_qc_input <- function(
    sample,
    meta
) {
  
  demux_rds <- file.path(
    outdir,
    "02_demuxed",
    paste0(
      sample,
      "_demuxed.rds"
    )
  )
  
  loaded_rds <- file.path(
    outdir,
    "01_loaded",
    paste0(
      sample,
      "_raw.rds"
    )
  )
  
  if (
    "demux" %in% steps &&
    sample_has_demux(meta)
  ) {
    return(demux_rds)
  }
  
  if (
    identical(steps, "qc") &&
    sample_has_demux(meta) &&
    file.exists(demux_rds)
  ) {
    
    msg(
      sample,
      ":: using existing demultiplexed object for QC"
    )
    
    return(demux_rds)
  }
  
  loaded_rds
}


# ============================================================ #
# Processing state
# ============================================================ #

successful_samples <- character()
failed_samples <- character()

per_sample_steps <- c(
  "load",
  "demux",
  "qc",
  "contam"
)


# ============================================================ #
# Per-sample processing
# ============================================================ #

for (sample in samples) {
  
  meta <- master[
    master$SampleID == sample,
    ,
    drop = FALSE
  ]
  
  sample_success <- TRUE
  
  sample_steps <- intersect(
    steps,
    per_sample_steps
  )
  
  for (step in sample_steps) {
    
    msg(
      "──",
      sample,
      "::",
      step
    )
    
    step_success <- tryCatch(
      expr = {
        
        # ---------------------------------------------------- #
        # Load Cell Ranger matrix
        # ---------------------------------------------------- #
        
        if (identical(step, "load")) {
          
          source_pipeline_script(
            "01_load_seurat.R"
          )
          
          cellranger_dir <- resolve_project_path(
            get_col(
              meta,
              "CellRanger_Dir"
            )
          )
          
          if (
            is.null(cellranger_dir) ||
            !dir.exists(cellranger_dir)
          ) {
            
            fail(
              "Cell Ranger directory not found:",
              cellranger_dir %||%
                "<not supplied>"
            )
          }
          
          # Confirm that a raw matrix exists for downstream
          # ambient-RNA decontamination.
          resolve_raw_matrix(
            cellranger_dir
          )
          
          run_load(
            sample = sample,
            cellranger = cellranger_dir,
            
            species = get_col(
              meta,
              "Species",
              cfg$species %||%
                "human"
            ),
            
            use_filtered = to_bool(
              cfg$use_filtered %||%
                TRUE,
              default = TRUE
            ),
            
            outdir = file.path(
              outdir,
              "01_loaded"
            )
          )
          
          # Apply fixed donor assignment for non-multiplexed
          # samples such as FN_PCW6 and FN_PCW8.
          donor_id <- get_col(
            meta,
            "DonorID",
            NA_character_
          )
          
          loaded_rds <- file.path(
            outdir,
            "01_loaded",
            paste0(
              sample,
              "_raw.rds"
            )
          )
          
          if (
            !is_missing_value(donor_id) &&
            file.exists(loaded_rds)
          ) {
            
            object <- readRDS(
              loaded_rds
            )
            
            object$status <- as.character(
              donor_id
            )
            
            object$donor_id <- as.character(
              donor_id
            )
            
            saveRDS(
              object,
              loaded_rds
            )
            
            msg(
              sample,
              ":: fixed donor assignment:",
              donor_id
            )
          }
        }
        
        
        # ---------------------------------------------------- #
        # Add Demuxafy or bulk-genotype assignments
        # ---------------------------------------------------- #
        
        if (identical(step, "demux")) {
          
          if (!sample_has_demux(meta)) {
            
            warn(
              sample,
              "has no Demuxafy or bulk-mapping input;",
              "demultiplexing skipped."
            )
            
          } else {
            
            source_pipeline_script(
              "02_add_demux.R"
            )
            
            loaded_rds <- file.path(
              outdir,
              "01_loaded",
              paste0(
                sample,
                "_raw.rds"
              )
            )
            
            if (!file.exists(loaded_rds)) {
              
              fail(
                "Loaded Seurat object not found:",
                loaded_rds,
                "\nRun the load step first."
              )
            }
            
            demux_inputs <- resolve_demux_inputs(
              sample = sample,
              meta = meta
            )
            
            demuxafy_dir <-
              demux_inputs$demuxafy_dir
            
            bulk_tsv <-
              demux_inputs$bulk_tsv
            
            if (
              is.null(demuxafy_dir) &&
              is.null(bulk_tsv)
            ) {
              
              fail(
                paste0(
                  "No valid Demuxafy or bulk-mapping ",
                  "input found for:"
                ),
                sample
              )
            }
            
            demux_cfg <- cfg$demux %||%
              list()
            
            assignment_preference <-
              demux_cfg$assignment_preference %||%
              cfg$demux_assignment_preference %||%
              "bulk"
            
            # Bulk assignment cannot be requested without a
            # valid bulk-mapping table.
            if (
              identical(
                assignment_preference,
                "bulk"
              ) &&
              is.null(bulk_tsv)
            ) {
              
              fail(
                paste0(
                  "assignment_preference='bulk', but no ",
                  "valid BulkMapping_TSV was found for:"
                ),
                sample
              )
            }
            
            run_demux(
              sample = sample,
              rds_in = loaded_rds,
              demuxafy_dir = demuxafy_dir,
              bulk_tsv = bulk_tsv,
              
              remove_doublets = to_bool(
                demux_cfg$remove_doublets %||%
                  cfg$remove_demux_doublets %||%
                  FALSE,
                default = FALSE
              ),
              
              remove_unassigned = to_bool(
                demux_cfg$remove_unassigned %||%
                  cfg$remove_demux_unassigned %||%
                  TRUE,
                default = TRUE
              ),
              
              assignment_preference =
                assignment_preference,
              
              keep_original_metadata = to_bool(
                demux_cfg$keep_original_metadata %||%
                  cfg$keep_demux_original_metadata %||%
                  FALSE,
                default = FALSE
              ),
              
              metadata_mode =
                demux_cfg$metadata_mode %||%
                cfg$demux_metadata_mode %||%
                "full",
              
              outdir = file.path(
                outdir,
                "02_demuxed"
              ),
              
              verbose = TRUE
            )
          }
        }
        
        
        # ---------------------------------------------------- #
        # Quality control
        # ---------------------------------------------------- #
        
        if (identical(step, "qc")) {
          
          source_pipeline_script(
            "03_run_qc.R"
          )
          
          qc_input <- get_qc_input(
            sample = sample,
            meta = meta
          )
          
          if (
            is.null(qc_input) ||
            !file.exists(qc_input)
          ) {
            
            fail(
              "QC input not found:",
              qc_input %||%
                "<not supplied>"
            )
          }
          
          run_qc_step(
            sample = sample,
            rds_in = qc_input,
            
            min_feat = as.integer(
              get_col(
                meta,
                "Min_Feat",
                cfg$default_min_feat %||%
                  200
              )
            ),
            
            min_umi = as.integer(
              get_col(
                meta,
                "Min_UMI",
                cfg$default_min_umi %||%
                  500
              )
            ),
            
            max_mito = as.numeric(
              get_col(
                meta,
                "Max_Mito",
                cfg$default_max_mito %||%
                  5
              )
            ),
            
            mad_n = as.integer(
              get_col(
                meta,
                "MAD_N",
                cfg$default_mad_n %||%
                  5
              )
            ),
            
            rm_dbl = to_bool(
              get_col(
                meta,
                "Rm_Dbl",
                cfg$default_rm_dbl %||%
                  FALSE
              ),
              default = FALSE
            ),
            
            dbl_score = as.numeric(
              get_col(
                meta,
                "Dbl_Score",
                cfg$default_dbl_score %||%
                  0.5
              )
            ),
            
            species = get_col(
              meta,
              "Species",
              cfg$species %||%
                "human"
            ),
            
            outdir = file.path(
              outdir,
              "03_qc"
            )
          )
        }
        
        
        # ---------------------------------------------------- #
        # Ambient RNA decontamination
        # ---------------------------------------------------- #
        
        if (identical(step, "contam")) {
          
          source_pipeline_script(
            "04_decontaminate.R"
          )
          
          qc_rds <- file.path(
            outdir,
            "03_qc",
            paste0(
              sample,
              "_qc.rds"
            )
          )
          
          if (!file.exists(qc_rds)) {
            
            fail(
              "QC object not found:",
              qc_rds,
              "\nRun the QC step first, for example:",
              paste(
                "Rscript run_pipeline.R",
                "--step qc",
                "--sample",
                sample
              )
            )
          }
          
          cellranger_dir <- resolve_project_path(
            get_col(
              meta,
              "CellRanger_Dir"
            )
          )
          
          if (
            is.null(cellranger_dir) ||
            !dir.exists(cellranger_dir)
          ) {
            
            fail(
              "Cell Ranger directory not found:",
              cellranger_dir %||%
                "<not supplied>"
            )
          }
          
          raw_input <- resolve_raw_matrix(
            cellranger_dir
          )
          
          run_decontam(
            sample = sample,
            rds_in = qc_rds,
            raw_dir = raw_input,
            
            assays = cfg$contam_assays %||%
              c(
                "RNA",
                "DecontX"
              ),
            
            preferred_assay =
              cfg$preferred_contam_assay %||%
              "DecontX",
            
            create_clean_alias = to_bool(
              cfg$create_clean_assay %||%
                TRUE,
              default = TRUE
            ),
            
            normalize_corrected = to_bool(
              cfg$normalize_corrected_assay %||%
                FALSE,
              default = FALSE
            ),
            
            outdir = file.path(
              outdir,
              "04_clean"
            )
          )
        }
        
        TRUE
      },
      
      error = function(error) {
        
        cat(
          sprintf(
            "[ERROR] %s :: %s — %s\n",
            sample,
            step,
            conditionMessage(error)
          )
        )
        
        FALSE
      }
    )
    
    if (!isTRUE(step_success)) {
      
      sample_success <- FALSE
      
      if (!isTRUE(
        options$continue_on_error
      )) {
        
        fail(
          "Pipeline stopped after failure:",
          sample,
          "::",
          step
        )
      }
      
      break
    }
    
    invisible(
      gc()
    )
  }
  
  if (isTRUE(sample_success)) {
    
    successful_samples <- c(
      successful_samples,
      sample
    )
    
  } else {
    
    failed_samples <- c(
      failed_samples,
      sample
    )
  }
}


# ============================================================ #
# Integration
# ============================================================ #

if ("integrate" %in% steps) {
  
  msg(
    "── integration"
  )
  
  source_pipeline_script(
    "05_integrate.R"
  )
  
  integration_samples <- if (
    length(
      intersect(
        steps,
        per_sample_steps
      )
    ) > 0L
  ) {
    
    successful_samples
    
  } else {
    
    samples
  }
  
  if (
    is.null(options$input_rds) &&
    length(integration_samples) == 0L
  ) {
    
    fail(
      "No successfully processed samples",
      "are available for integration."
    )
  }
  
  rds_list <- if (!is.null(
    options$input_rds
  )) {
    
    split_paths(
      options$input_rds
    )
    
  } else if ("contam" %in% steps) {
    
    file.path(
      outdir,
      "04_clean",
      paste0(
        integration_samples,
        "_clean.rds"
      )
    )
    
  } else {
    
    file.path(
      outdir,
      "03_qc",
      paste0(
        integration_samples,
        "_qc.rds"
      )
    )
  }
  
  rds_list <- vapply(
    rds_list,
    resolve_project_path,
    character(1L)
  )
  
  missing_inputs <- rds_list[
    !file.exists(rds_list)
  ]
  
  if (length(missing_inputs) > 0L) {
    
    fail(
      "Missing integration input(s):",
      paste(
        missing_inputs,
        collapse = ", "
      )
    )
  }
  
  run_integration(
    rds_paths = rds_list,
    method = cfg$integration_method %||% "Harmony", 
    batch_col = cfg$integration_batch_col %||% "status",
    dims = as.integer(cfg$pca_dims %||% 30),
    resolution = as.numeric(cfg$cluster_resolution %||% 0.5),
    nfeatures = as.integer(cfg$nfeatures %||% 3000),
    npcs = as.integer(cfg$npcs %||% 50),
    seed = as.integer(cfg$seed %||% 1234),
    outdir = file.path(outdir, "05_integrated")
  )
  
  invisible(
    gc()
  )
}


# ============================================================ #
# Differential expression
# ============================================================ #

if ("deg" %in% steps) {
  
  msg(
    "── DEG analysis"
  )
  
  integrated_rds <- if (!is.null(
    options$input_rds
  )) {
    
    supplied <- split_paths(
      options$input_rds
    )
    
    if (length(supplied) != 1L) {
      
      fail(
        "--input_rds must contain exactly one",
        "integrated RDS for the DEG step."
      )
    }
    
    resolve_project_path(
      supplied[[1L]]
    )
    
  } else {
    
    file.path(
      outdir,
      "05_integrated",
      "integrated.rds"
    )
  }
  
  if (
    is.null(integrated_rds) ||
    !file.exists(integrated_rds)
  ) {
    
    fail(
      "Integrated RDS not found:",
      integrated_rds %||%
        "<not supplied>"
    )
  }
  
  contrasts_file <- resolve_required_path(
    cfg$contrasts_file,
    "Contrasts file"
  )
  
  deg_outdir <- file.path(
    outdir,
    "06_deg"
  )
  
  dir.create(
    deg_outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  source_pipeline_script(
    "06_run_deg.R"
  )
  
  if (!exists(
    "run_deg",
    mode = "function"
  )) {
    
    fail(
      "06_run_deg.R does not define",
      "a run_deg() function."
    )
  }
  
  deg_results <- run_deg(
    rds_in = integrated_rds,
    contrasts_file = contrasts_file,
    
    method = cfg$deg_method %||%
      "pseudo_bulk",
    
    donor_col = cfg$donor_col %||%
      "donor_id",
    
    fdr_thr = as.numeric(
      cfg$fdr_threshold %||%
        0.05
    ),
    
    lfc_thr = as.numeric(
      cfg$lfc_threshold %||%
        0.5
    ),
    
    top_n = as.integer(
      cfg$top_n_genes %||%
        20
    ),
    
    run_pathway = to_bool(
      cfg$run_pathway %||%
        TRUE,
      default = TRUE
    ),
    
    species = cfg$species %||%
      "human",
    
    outdir = deg_outdir
  )
  
  msg(
    "DEG analysis complete:",
    deg_outdir
  )
  
  invisible(
    gc()
  )
}


# ============================================================ #
# Final status
# ============================================================ #

if (length(successful_samples) > 0L) {
  
  msg(
    "Successful sample(s):",
    paste(
      unique(successful_samples),
      collapse = ", "
    )
  )
}

if (length(failed_samples) > 0L) {
  
  warn(
    "Failed sample(s):",
    paste(
      unique(failed_samples),
      collapse = ", "
    )
  )
}

msg(
  "Pipeline run complete."
)