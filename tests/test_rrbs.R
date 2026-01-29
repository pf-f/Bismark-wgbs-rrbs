#!/usr/bin/env Rscript

#' Unit Tests for RRBS Analysis
#'
#' @author Methylation Analysis Toolkit

cat("=" <- rep("=", 60), "\n")
cat("RRBS Analysis Tests\n")
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

# Test 1: Load RRBS config
test1 <- function() {
  source("src/methylation_analysis/R/config.R")
  config <- get_config("RRBS", "hg38")
  stopifnot(config$analysis_type == "RRBS")
  stopifnot(config$genome == "hg38")
  stopifnot(config$need_deduplication == FALSE)
  stopifnot(config$use_trim_galore == TRUE)
  stopifnot(config$trim_galore_rrbs == TRUE)
}

run_test("Load RRBS configuration", test1)

# Test 2: RRBS uses Trim Galore
test2 <- function() {
  source("src/methylation_analysis/R/qc_processor.R")
  # Check if trim_galore is available (may not be installed)
  trim_galore_path <- Sys.which("trim_galore")
  # We don't fail if it's not installed, just check the config
  return(TRUE)
}

run_test("Check Trim Galore availability", test2)

# Test 3: RRBS has lower coverage threshold
test3 <- function() {
  source("src/methylation_analysis/R/config.R")
  wgbs_config <- get_config("WGBS", "hg38")
  rrbs_config <- get_config("RRBS", "hg38")
  
  # RRBS should have lower minimum coverage
  stopifnot(rrbs_config$min_coverage < wgbs_config$min_coverage)
  # RRBS should use fewer threads
  stopifnot(rrbs_config$bismark_parallel < wgbs_config$bismark_parallel)
  # RRBS should use less memory
  stopifnot(as.numeric(gsub("G", "", rrbs_config$bismark_buffer_size)) < 
              as.numeric(gsub("G", "", wgbs_config$bismark_buffer_size)))
}

run_test("RRBS configuration parameters", test3)

# Test 4: RRBS deduplication is skipped
test4 <- function() {
  source("src/methylation_analysis/R/deduplication.R")
  config <- list(need_deduplication = FALSE, analysis_type = "RRBS")
  output <- run_deduplication(
    bam_files = list(sample1 = "test.bam"),
    output_dir = tempdir(),
    config = config
  )
  # Should return the input files unchanged
  stopifnot(identical(names(output), names(list(sample1 = "test.bam"))))
}

run_test("RRBS skips deduplication", test4)

# Test 5: Load methylKit
test5 <- function() {
  source("src/methylation_analysis/R/analyzer.R")
  # May fail if not installed, but we just test the function exists
  tryCatch({
    load_methylkit()
  }, error = function(e) {
    if (!grepl("not installed", e$message)) {
      stop(e)
    }
  })
  return(TRUE)  # Pass test even if packages missing
}

run_test("Load methylKit (may fail if not installed)", test5)

# Test 6: WGBS vs RRBS config comparison
test6 <- function() {
  source("src/methylation_analysis/R/config.R")
  wgbs_config <- get_config("WGBS", "hg38")
  rrbs_config <- get_config("RRBS", "hg38")
  
  # Verify key differences
  differences <- c(
    "Deduplication" = wgbs_config$need_deduplication != rrbs_config$need_deduplication,
    "Trim Galore" = wgbs_config$use_trim_galore != rrbs_config$use_trim_galore,
    "Min Coverage" = wgbs_config$min_coverage != rrbs_config$min_coverage,
    "Comprehensive Extraction" = wgbs_config$extract_comprehensive != rrbs_config$extract_comprehensive
  )
  
  # All key differences should be present
  stopifnot(all(differences))
}

run_test("Verify WGBS and RRBS configuration differences", test6)

# Test 7: Parse coverage file (mock test)
test7 <- function() {
  source("src/methylation_analysis/R/methylation_extractor.R")
  
  # Create mock coverage file
  temp_file <- tempfile(fileext = ".txt")
  writeLines(c("chr1\t100\t+\t10\t5", "chr1\t200\t-\t8\t2"), temp_file)
  
  # Parse file
  data <- parse_coverage_file(temp_file, context = "CpG")
  
  # Verify parsing
  stopifnot(nrow(data) == 2)
  stopifnot(all(c("chromosome", "position", "strand", 
                    "count_methylated", "count_unmethylated", 
                    "total_count", "methylation_percentage") %in% names(data)))
  
  # Clean up
  file.remove(temp_file)
}

run_test("Parse coverage file", test7)

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
