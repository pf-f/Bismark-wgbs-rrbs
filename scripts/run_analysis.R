#!/usr/bin/env Rscript

#' Methylation Analysis Pipeline - Command Line Entry Point
#'
#' This script provides a command-line interface for running the complete
#' methylation analysis pipeline for both WGBS and RRBS data.
#'
#' @author Methylation Analysis Toolkit
#' @version 1.0.0

suppressPackageStartupMessages({
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-t", "--type"), 
              type = "character", 
              default = "WGBS",
              help = "Analysis type: WGBS or RRBS [default: %default]"),
  
  make_option(c("-g", "--genome"), 
              type = "character", 
              default = "hg38",
              help = "Genome version (e.g., hg38, mm10) [default: %default]"),
  
  make_option(c("-i", "--input"), 
              type = "character",
              default = "data/raw_fastq",
              help = "Input directory with FASTQ files [default: %default]"),
  
  make_option(c("-c", "--config"), 
              type = "character", 
              default = NULL,
              help = "Configuration file (YAML). If NULL, uses default config"),
  
  make_option(c("-o", "--output"), 
              type = "character", 
              default = "outputs",
              help = "Output directory [default: %default]"),
  
  make_option(c("-s", "--steps"), 
              type = "character", 
              default = NULL,
              help = "Comma-separated steps to run: qc,align,dedup,extract,diff,annotate,viz"),
  
  make_option(c("--threads"), 
              type = "integer", 
              default = 8,
              help = "Number of threads for parallel processing [default: %default]"),
  
  make_option(c("--sample-info"), 
              type = "character", 
              default = NULL,
              help = "Sample information file (CSV with columns: run_accession, sample_alias, treatment)"),
  
  make_option(c("--min-cov"), 
              type = "integer", 
              default = 10,
              help = "Minimum coverage threshold [default: %default]"),
  
  make_option(c("--diff-thresh"), 
              type = "integer", 
              default = 25,
              help = "Methylation difference threshold for DMRs (%) [default: %default]"),
  
  make_option(c("--qvalue"), 
              type = "double", 
              default = 0.01,
              help = "Q-value threshold for significance [default: %default]"),
  
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = TRUE,
              help = "Print verbose output [default: TRUE]"),
  
  make_option(c("--version"), 
              action = "store_true", 
              default = FALSE,
              help = "Print version and exit")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "\nMethylation Analysis Pipeline for WGBS and RRBS Data")
opt <- parse_args(opt_parser)

# Check for version flag
if (opt$version) {
  cat("Methylation Analysis Pipeline v1.0.0\n")
  cat("Supports WGBS and RRBS analysis workflows\n")
  quit(status = 0)
}

# Validate analysis type
if (!opt$type %in% c("WGBS", "RRBS")) {
  print_help(opt_parser)
  cat("\nError: Invalid analysis type. Must be 'WGBS' or 'RRBS'.\n")
  quit(status = 1)
}

# Print header
cat("\n")
cat("=" <- rep("=", 70), "\n")
cat("Methylation Analysis Pipeline\n")
cat("=" <- rep("=", 70), "\n")
cat("Analysis Type:", opt$type, "\n")
cat("Genome:", opt$genome, "\n")
cat("Input Directory:", opt$input, "\n")
cat("Output Directory:", opt$output, "\n")
cat("=" <- rep("=", 70), "\n\n")

# Load modules
module_dir <- file.path(dirname(getwd()), "src", "methylation_analysis", "R")
if (!dir.exists(module_dir)) {
  # Try package installation path
  module_dir <- system.file("R", package = "methylation_analysis", mustWork = FALSE)
  if (module_dir == "") {
    stop("Cannot find R modules directory: ", module_dir)
  }
}

# Source all modules
module_files <- list.files(module_dir, pattern = "\\.R$", full.names = TRUE)
for (module_file in module_files) {
  source(module_file)
}

# Load configuration
if (!is.null(opt$config)) {
  cat("Loading configuration from:", opt$config, "\n")
  config <- tryCatch({
    load_config(opt$config)
  }, error = function(e) {
    stop("Error loading configuration: ", e$message)
  })
} else {
  cat("Using default configuration\n")
  config <- get_config(analysis_type = opt$type, genome = opt$genome)
}

# Override config with command line options
config$output_dir <- opt$output
config$threads <- opt$threads
config$min_coverage <- opt$min_cov
config$diff_meth_threshold <- opt$diff_thresh
config$diff_qvalue_threshold <- opt$qvalue

# Update output directories
config$output_dirs <- list(
  clip_reads = file.path(opt$output, "01_clip_reads"),
  meth_bam = file.path(opt$output, "02_meth_bam"),
  dedu_meth = file.path(opt$output, "03_dedu_meth"),
  methinfo = file.path(opt$output, "04_methinfo"),
  results = file.path(opt$output, "05_results")
)

# Load sample information if provided
sample_info <- NULL
if (!is.null(opt$sample_info)) {
  if (!file.exists(opt$sample_info)) {
    stop("Sample info file not found: ", opt$sample_info)
  }
  
  cat("Loading sample information from:", opt$sample_info, "\n")
  sample_info <- read.csv(opt$sample_info, stringsAsFactors = FALSE)
  
  # Add treatment code if not present
  if (!"treatment_code" %in% names(sample_info)) {
    # Map treatment to code
    treatment_map <- c("CTC" = 2, "T" = 1, "Normal" = 0)
    sample_info$treatment_code <- treatment_map[as.character(sample_info$treatment)]
  }
  
  cat("  Loaded", nrow(sample_info), "samples\n")
}

# Parse steps
if (!is.null(opt$steps)) {
  steps <- strsplit(opt$steps, ",")[[1]]
  steps <- trimws(steps)
  steps <- toupper(steps)
} else {
  steps <- NULL
}

# Create log directory
log_dir <- file.path(opt$output, "log")
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

# Set up logging
if (opt$verbose) {
  log_file <- file.path(log_dir, paste0("pipeline_", Sys.Date(), ".log"))
  cat("Log file:", log_file, "\n\n")
  
  # Redirect output to log file
  sink(log_file, split = TRUE)
  
  # Print configuration
  print_config(config)
  
  sink()
}

# Run pipeline
cat("\nStarting analysis pipeline...\n\n")

start_time <- Sys.time()

tryCatch({
  results <- run_pipeline(
    config = config,
    input_fastq_dir = opt$input,
    sample_info = sample_info,
    steps = steps
  )
  
  # Generate summary report
  generate_pipeline_summary(results, config)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat("\n")
  cat("=" <- rep("=", 70), "\n")
  cat("Pipeline completed successfully!\n")
  cat("Total time:", round(as.numeric(elapsed), 2), "minutes\n")
  cat("=" <- rep("=", 70), "\n")
  cat("\nOutput directory:", opt$output, "\n")
  cat("Log directory:", log_dir, "\n")
  
}, error = function(e) {
  cat("\n")
  cat("=" <- rep("=", 70), "\n")
  cat("ERROR: Pipeline failed!\n")
  cat("=" <- rep("=", 70), "\n")
  cat("Message:", e$message, "\n\n")
  
  # Print call stack
  if (opt$verbose) {
    cat("Call stack:\n")
    print(e$call)
  }
  
  quit(status = 1)
})

# Quit successfully
quit(status = 0)
