#' Configuration Management for Methylation Analysis
#'
#' This module handles configuration for both WGBS and RRBS workflows
#'
#' @author Methylation Analysis Toolkit
#' @export

# Constants
VALID_ANALYSIS_TYPES <- c("WGBS", "RRBS")

#' Get default configuration for analysis type
#'
#' @param analysis_type Character: "WGBS" or "RRBS"
#' @param genome Character: genome version (e.g., "hg38", "mm10")
#' @return List containing configuration parameters
#' @export
get_config <- function(analysis_type = "WGBS", genome = "hg38") {
  
  # Validate analysis type
  if (!analysis_type %in% VALID_ANALYSIS_TYPES) {
    stop("Invalid analysis_type. Must be one of: ", 
         paste(VALID_ANALYSIS_TYPES, collapse = ", "))
  }
  
  # Common configuration
  common_config <- list(
    # Directories
    ref_dir = "data/ref",
    output_dir = "outputs",
    log_dir = "log",
    
    # Genome
    genome = genome,
    
    # Processing threads
    threads = 8,
    
    # File patterns
    fastq_pattern = "*_\\d.fastq.gz",
    
    # Quality control
    min_read_length = 75,
    trim_sliding_window = "4:15",
    
    # Coverage thresholds
    min_coverage = 10,
    max_coverage_percentile = 99.9,
    
    # Methylation analysis
    read_context = "CpG",
    diff_meth_threshold = 25,
    diff_qvalue_threshold = 0.01
  )
  
  # WGBS-specific configuration
  wgbs_config <- list(
    # Deduplication
    need_deduplication = TRUE,
    
    # Bismark parameters
    bismark_parallel = 8,
    bismark_buffer_size = "20G",
    
    # Methylation extraction
    extract_comprehensive = TRUE,
    extract_bedgraph = TRUE,
    extract_cytosine_report = TRUE
  )
  
  # RRBS-specific configuration
  rrbs_config <- list(
    # Deduplication (skip for RRBS)
    need_deduplication = FALSE,
    
    # Use trim-galore for RRBS
    use_trim_galore = TRUE,
    
    # RRBS-specific parameters
    trim_galore_rrbs = TRUE,
    
    # Lower coverage threshold for RRBS
    min_coverage = 5,
    
    # Bismark parameters
    bismark_parallel = 4,
    bismark_buffer_size = "10G",
    
    # Methylation extraction
    extract_comprehensive = FALSE,
    extract_bedgraph = TRUE,
    extract_cytosine_report = TRUE
  )
  
  # Merge configurations
  config <- if (analysis_type == "WGBS") {
    c(common_config, wgbs_config)
  } else {
    c(common_config, rrbs_config)
  }
  
  # Add analysis type
  config$analysis_type <- analysis_type
  
  # Add output subdirectories
  config$output_dirs <- list(
    clip_reads = file.path(config$output_dir, "01_clip_reads"),
    meth_bam = file.path(config$output_dir, "02_meth_bam"),
    dedu_meth = file.path(config$output_dir, "03_dedu_meth"),
    methinfo = file.path(config$output_dir, "04_methinfo"),
    results = file.path(config$output_dir, "05_results")
  )
  
  return(config)
}

#' Load configuration from YAML file
#'
#' @param config_file Path to configuration file
#' @return List containing configuration parameters
#' @export
load_config <- function(config_file) {
  
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }
  
  # Check if yaml package is available
  if (!requireNamespace("yaml", quietly = TRUE)) {
    # Fallback: read as text and parse manually
    message("yaml package not available, using basic configuration")
    return(get_config())
  }
  
  config <- yaml::read_yaml(config_file)
  
  # Validate required fields
  required_fields <- c("analysis_type", "genome")
  missing_fields <- setdiff(required_fields, names(config))
  
  if (length(missing_fields) > 0) {
    stop("Missing required fields in config: ", 
         paste(missing_fields, collapse = ", "))
  }
  
  # Merge with default configuration
  default_config <- get_config(config$analysis_type, config$genome)
  config <- modifyList(default_config, config)
  
  return(config)
}

#' Print configuration summary
#'
#' @param config Configuration list
#' @export
print_config <- function(config) {
  
  cat("=" <- rep("=", 60), "\n")
  cat("Methylation Analysis Configuration\n")
  cat("=" <- rep("=", 60), "\n\n")
  
  cat("Analysis Type:", config$analysis_type, "\n")
  cat("Genome:", config$genome, "\n")
  cat("Threads:", config$threads, "\n\n")
  
  cat("Parameters:\n")
  cat("  Min Coverage:", config$min_coverage, "\n")
  cat("  Min Read Length:", config$min_read_length, "\n")
  cat("  Diff Meth Threshold:", config$diff_meth_threshold, "%\n")
  cat("  Diff Q-value Threshold:", config$diff_qvalue_threshold, "\n")
  
  cat("\nDeduplication:", 
      ifelse(config$need_deduplication, "Yes", "No"), "\n")
  
  cat("\nOutput Directories:\n")
  for (dir_name in names(config$output_dirs)) {
    cat("  ", dir_name, ":", config$output_dirs[[dir_name]], "\n")
  }
  
  cat("\n", "=" <- rep("=", 60), "\n")
  
  invisible(config)
}

#' Validate configuration
#'
#' @param config Configuration list
#' @return TRUE if valid, otherwise stops with error
#' @export
validate_config <- function(config) {
  
  # Check required directories exist
  if (!dir.exists(config$ref_dir)) {
    stop("Reference directory not found: ", config$ref_dir)
  }
  
  # Check genome exists
  genome_path <- file.path(config$ref_dir, toupper(config$genome))
  if (!dir.exists(genome_path)) {
    stop("Genome directory not found: ", genome_path)
  }
  
  # Create output directories
  for (dir_name in names(config$output_dirs)) {
    dir_path <- config$output_dirs[[dir_name]]
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Created directory: ", dir_path)
    }
  }
  
  return(TRUE)
}
