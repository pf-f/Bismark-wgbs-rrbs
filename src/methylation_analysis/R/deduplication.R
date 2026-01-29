#' Deduplication Module for WGBS
#'
#' This module handles read deduplication for WGBS data
#' Note: RRBS data should NOT be deduplicated
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Run Bismark deduplication
#'
#' @param bam_files Character vector of BAM files
#' @param output_dir Output directory for deduplicated BAMs
#' @param config Configuration list
#' @return List of deduplicated BAM file paths
#' @export
run_deduplication <- function(bam_files, output_dir, config) {
  
  # Check if deduplication is needed
  if (!config$need_deduplication) {
    message("Deduplication skipped for ", config$analysis_type, " analysis.")
    message("Deduplication is only required for WGBS, not RRBS.")
    return(bam_files)
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if deduplicate_bismark is available
  dedup_path <- Sys.which("deduplicate_bismark")
  if (dedup_path == "") {
    stop("deduplicate_bismark not found. Please ensure bismark is properly installed.")
  }
  
  message("=" <- rep("=", 60))
  message("Running Bismark Deduplication")
  message("=" <- rep("=", 60))
  message("Processing ", length(bam_files), " BAM files...")
  
  dedup_bams <- list()
  
  for (bam_file in bam_files) {
    sample_name <- gsub("\\.bam$", "", basename(bam_file))
    sample_name <- gsub("_bismark.*", "", sample_name)
    
    message("\nDeduplicating sample: ", sample_name)
    
    # Build command
    cmd <- paste(
      dedup_path,
      "--bam",
      "-p",
      "-o", sample_name, "_pe_dedup.bam",
      "--output_dir", output_dir,
      bam_file
    )
    
    # Execute command
    start_time <- Sys.time()
    result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    end_time <- Sys.time()
    
    if (result != 0) {
      warning("Deduplication failed for sample: ", sample_name)
    } else {
      elapsed <- difftime(end_time, start_time, units = "mins")
      message("  Completed in", round(as.numeric(elapsed), 2), "minutes")
      
      # Find output file
      dedup_pattern <- paste0(sample_name, "_pe_dedup.*\\.bam$")
      dedup_files <- list.files(output_dir, pattern = dedup_pattern, full.names = TRUE)
      
      if (length(dedup_files) > 0) {
        dedup_bams[[sample_name]] <- dedup_files[1]
        message("  Output BAM:", dedup_files[1])
      }
    }
  }
  
  message("\nDeduplication completed!")
  message("Successfully deduplicated:", length(dedup_bams), "/", length(bam_files), "samples")
  
  return(dedup_bams)
}

#' Check deduplication report
#'
#' @param bam_file Path to BAM file
#' @return List with deduplication statistics
#' @export
check_dedup_stats <- function(bam_file) {
  
  # Check if samtools is available
  samtools_path <- Sys.which("samtools")
  if (samtools_path == "") {
    stop("samtools not found. Please install: conda install -c bioconda samtools")
  }
  
  # Count total reads
  total_reads <- as.numeric(system2(samtools_path, 
                                    c("view", "-c", bam_file), 
                                    stdout = TRUE))
  
  # Count properly paired reads
  paired_reads <- as.numeric(system2(samtools_path, 
                                      c("view", "-c", "-f", "2", bam_file), 
                                      stdout = TRUE))
  
  # Count duplicate reads (flag 1024)
  duplicate_reads <- as.numeric(system2(samtools_path, 
                                        c("view", "-c", "-f", "1024", bam_file), 
                                        stdout = TRUE))
  
  stats <- list(
    total_reads = total_reads,
    paired_reads = paired_reads,
    duplicate_reads = duplicate_reads,
    duplication_rate = duplicate_reads / total_reads * 100
  )
  
  return(stats)
}

#' Print deduplication statistics
#'
#' @param stats_list List of statistics for multiple samples
#' @export
print_dedup_stats <- function(stats_list) {
  
  if (length(stats_list) == 0) {
    message("No deduplication statistics available.")
    return(invisible(NULL))
  }
  
  # Create data frame
  df <- do.call(rbind, lapply(names(stats_list), function(sample_name) {
    stats <- stats_list[[sample_name]]
    data.frame(
      Sample = sample_name,
      Total_Reads = stats$total_reads,
      Paired_Reads = stats$paired_reads,
      Duplicate_Reads = stats$duplicate_reads,
      Duplication_Rate = paste0(round(stats$duplication_rate, 2), "%")
    )
  }))
  
  cat("\nDeduplication Statistics:\n")
  print(df, row.names = FALSE)
  cat("\n")
  
  invisible(df)
}

#' Remove temporary Bismark files
#'
#' @param bismark_dir Directory containing Bismark output
#' @export
cleanup_bismark_temp <- function(bismark_dir) {
  
  # Remove temporary alignment files
  temp_files <- list.files(bismark_dir, pattern = "\\.temp\\.")
  
  if (length(temp_files) > 0) {
    message("Removing ", length(temp_files), " temporary files...")
    file.remove(file.path(bismark_dir, temp_files))
  }
  
  # Remove intermediate files
  intermediate_patterns <- c("*.bt2$", "*\\.fa$", "*\\_val_\\d+\\.fq\\.gz$")
  
  for (pattern in intermediate_patterns) {
    files <- list.files(bismark_dir, pattern = pattern, full.names = TRUE)
    if (length(files) > 0) {
      message("Removing ", length(files), " intermediate files matching: ", pattern)
      file.remove(files)
    }
  }
  
  message("Cleanup completed.")
  
  invisible(NULL)
}

#' Estimate required memory for alignment
#'
#' @param genome_size Genome size in bases
#' @param threads Number of threads
#' @return Estimated memory in GB
#' @export
estimate_alignment_memory <- function(genome_size = 3e9, threads = 4) {
  
  # Bismark memory estimation
  # Based on Bismark documentation
  
  bowtie2_memory <- 4  # GB (small index)
  bismark_overhead <- 1  # GB
  
  # Additional memory for parallel processing
  parallel_memory <- 0.5 * threads
  
  total_memory <- bowtie2_memory + bismark_overhead + parallel_memory
  
  message("Estimated memory requirement:", round(total_memory, 2), "GB")
  message("  Bowtie2 index:", bowtie2_memory, "GB")
  message("  Bismark overhead:", bismark_overhead, "GB")
  message("  Parallel processing:", parallel_memory, "GB")
  
  return(total_memory)
}

#' Validate deduplicated BAM files
#'
#' @param bam_files Character vector of BAM file paths
#' @return Logical vector indicating valid files
#' @export
validate_bam_files <- function(bam_files) {
  
  samtools_path <- Sys.which("samtools")
  if (samtools_path == "") {
    warning("samtools not found. Skipping BAM validation.")
    return(rep(TRUE, length(bam_files)))
  }
  
  valid <- sapply(bam_files, function(bam_file) {
    result <- system2(samtools_path, 
                       c("quickcheck", bam_file), 
                       stdout = FALSE, stderr = FALSE)
    return(result == 0)
  })
  
  names(valid) <- basename(bam_files)
  
  if (any(!valid)) {
    warning("Some BAM files failed validation:",
            paste(names(valid)[!valid], collapse = ", "))
  }
  
  return(valid)
}
