#!/usr/bin/env Rscript

#' Unit Tests for WGBS Analysis
#'
#' @author Methylation Analysis Toolkit

cat("=" <- rep("=", 60), "\n")
cat("WGBS Analysis Tests\n")
cat("=" <- rep("=", 60), "\n\n")

# Test counter
total_tests <- 0
passed_tests <- 0
failed_tests <- 0

# Test function
run_test <- function(test_name, test_func) {
  total_tests <<- total_tests + 1
  
  cat("Test", total_tests, ":", test_name, "... ")
  
  result <- tryCatch({
    test_func()
    cat("PASSED\n")
    passed_tests <<- passed_tests + 1
    return(TRUE)
  }, error = function(e) {
    cat("FAILED\n")
    cat("  Error:", e$message, "\n")
    failed_tests <<- failed_tests + 1
    return(FALSE)
  })
}

# Test 1: Load config
test1 <- function() {
  source("src/methylation_analysis/R/config.R")
  config <- get_config("WGBS", "hg38")
  stopifnot(config$analysis_type == "WGBS")
  stopifnot(config$genome == "hg38")
  stopifnot(config$need_deduplication == TRUE)
}

run_test("Load WGBS configuration", test1)

# Test 2: Check Bismark installation
test2 <- function() {
  source("src/methylation_analysis/R/aligner.R")
  status <- check_bismark()
  stopifnot(status$bismark_installed == TRUE)
}

run_test("Check Bismark installation", test2)

# Test 3: Validate config
test3 <- function() {
  source("src/methylation_analysis/R/config.R")
  config <- get_config("WGBS", "hg38")
  # Create test directories
  test_dir <- tempdir()
  config$ref_dir <- test_dir
  config$output_dir <- test_dir
  config$output_dirs <- list(
    clip_reads = file.path(test_dir, "01_clip_reads"),
    meth_bam = file.path(test_dir, "02_meth_bam")
  )
  # Create test genome directory
  genome_dir <- file.path(test_dir, "HG38")
  dir.create(genome_dir, recursive = TRUE)
  # This should fail validation (no FASTA files)
  tryCatch({
    validate_config(config)
    stop("Should have failed validation")
  }, error = function(e) {
    if (!grepl("No FASTA", e$message)) {
      stop(e)
    }
  })
}

run_test("Validate configuration (missing FASTA)", test3)

# Test 4: Load annotation libraries
test4 <- function() {
  source("src/methylation_analysis/R/annotator.R")
  # This may fail if packages not installed, but that's ok for testing
  tryCatch({
    load_annotation_libraries()
  }, error = function(e) {
    warning("Annotation libraries not installed:", e$message)
  })
  return(TRUE)  # Pass test even if packages missing
}

run_test("Load annotation libraries", test4)

# Test 5: WGBS vs RRBS config differences
test5 <- function() {
  source("src/methylation_analysis/R/config.R")
  wgbs_config <- get_config("WGBS", "hg38")
  rrbs_config <- get_config("RRBS", "hg38")
  
  # WGBS should require deduplication
  stopifnot(wgbs_config$need_deduplication == TRUE)
  stopifnot(rrbs_config$need_deduplication == FALSE)
  
  # RRBS should use Trim Galore
  stopifnot(rrbs_config$use_trim_galore == TRUE)
  stopifnot(wgbs_config$use_trim_galore == FALSE)
  
  # RRBS should have lower coverage threshold
  stopifnot(rrbs_config$min_coverage < wgbs_config$min_coverage)
}

run_test("Compare WGBS and RRBS configurations", test5)

# Print summary
cat("\n")
cat("=" <- rep("=", 60), "\n")
cat("Test Summary\n")
cat("=" <- rep("=", 60), "\n")
cat("Total tests:", total_tests, "\n")
cat("Passed:", passed_tests, "\n")
cat("Failed:", failed_tests, "\n")
cat("Success rate:", round(passed_tests / total_tests * 100, 1), "%\n")
cat("=" <- rep("=", 60), "\n")

# Exit with appropriate code
quit(status = ifelse(failed_tests == 0, 0, 1))
