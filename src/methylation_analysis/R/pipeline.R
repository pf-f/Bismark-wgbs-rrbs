#' Main Pipeline Controller
#'
#' This module orchestrates the complete methylation analysis workflow
#' for both WGBS and RRBS data
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Load all analysis modules
#'
#' @export
load_modules <- function() {
  
  source_dir <- system.file("R", package = "methylation_analysis", 
                          mustWork = FALSE)
  
  if (source_dir == "") {
    # Development mode: use relative path
    source_dir <- file.path(dirname(getwd()), "src", "methylation_analysis", "R")
  }
  
  # Source all R modules
  module_files <- list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
  
  for (module_file in module_files) {
    source(module_file)
  }
  
  message("Loaded", length(module_files), "analysis modules.")
  
  invisible(TRUE)
}

#' Main methylation analysis pipeline
#'
#' @param config Configuration list (from get_config or load_config)
#' @param input_fastq_dir Directory containing input FASTQ files
#' @param sample_info Data frame with sample information
#' @param steps Character vector of steps to run. If NULL, run all steps.
#'   Available steps: "qc", "align", "dedup", "extract", "diff", "annotate", "viz"
#' @return List with pipeline results
#' @export
run_pipeline <- function(config, input_fastq_dir, sample_info = NULL, 
                       steps = NULL) {
  
  # Load modules
  load_modules()
  
  # Validate configuration
  validate_config(config)
  
  # Print configuration
  print_config(config)
  
  # Define all available steps
  all_steps <- c("qc", "align", "dedup", "extract", "diff", "annotate", "viz")
  
  # If steps not specified, run all steps
  if (is.null(steps)) {
    steps <- all_steps
    message("\nRunning complete pipeline...")
  } else {
    # Validate requested steps
    invalid_steps <- setdiff(steps, all_steps)
    if (length(invalid_steps) > 0) {
      stop("Invalid steps: ", paste(invalid_steps, collapse = ", "),
           "\nAvailable steps: ", paste(all_steps, collapse = ", "))
    }
    message("\nRunning steps:", paste(steps, collapse = ", "))
  }
  
  # Initialize results list
  pipeline_results <- list()
  
  # Step 1: Quality Control
  if ("qc" %in% steps) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 1: Quality Control and Preprocessing")
    message("=" <- rep("=", 60))
    
    input_files <- list.files(input_fastq_dir, pattern = "\\.fastq\\.gz$", 
                            full.names = TRUE)
    
    if (length(input_files) == 0) {
      stop("No FASTQ files found in: ", input_fastq_dir)
    }
    
    pipeline_results$trim_output <- run_qc(
      input_files,
      output_dir = config$output_dir,
      config = config
    )
  }
  
  # Step 2: Alignment
  if ("align" %in% steps) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 2: Alignment with Bismark")
    message("=" <- rep("=", 60))
    
    # Prepare genome
    genome_path <- file.path(config$ref_dir, toupper(config$genome))
    bismark_genome_dir <- prepare_bismark_genome(
      genome_path,
      verbose = TRUE
    )
    
    # Get trimmed files
    if (is.null(pipeline_results$trim_output)) {
      # Try to find trimmed files in default location
      clip_dir <- config$output_dirs$clip_reads
      trimmed_files <- list.files(clip_dir, pattern = "_val_\\d\\.fq\\.gz|_clip\\.\\d\\.fq\\.gz$", 
                               full.names = TRUE)
      
      if (length(trimmed_files) == 0) {
        stop("No trimmed files found. Please run QC step first.")
      }
      
      # Group into pairs
      sample_names <- unique(gsub("_val_\\d\\.fq\\.gz|_clip\\.\\d\\.fq\\.gz$", "", 
                                 basename(trimmed_files)))
      trim_output <- lapply(sample_names, function(s) {
        r1 <- grep(paste0(s, "_val_1\\.fq\\.gz"), trimmed_files, value = TRUE)
        r2 <- grep(paste0(s, "_val_2\\.fq\\.gz"), trimmed_files, value = TRUE)
        list(r1 = ifelse(length(r1) > 0, r1, NA),
              r2 = ifelse(length(r2) > 0, r2, NA))
      })
      names(trim_output) <- sample_names
      pipeline_results$trim_output <- trim_output
    }
    
    # Run alignment
    pipeline_results$align_output <- run_bismark(
      pipeline_results$trim_output,
      genome_dir = bismark_genome_dir,
      output_dir = config$output_dirs$meth_bam,
      config = config
    )
  }
  
  # Step 3: Deduplication (WGBS only)
  if ("dedup" %in% steps && config$need_deduplication) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 3: Deduplication")
    message("=" <- rep("=", 60))
    
    pipeline_results$dedup_output <- run_deduplication(
      pipeline_results$align_output,
      output_dir = config$output_dirs$dedu_meth,
      config = config
    )
    
    # Sort BAM files for methylKit
    sorted_bams <- sort_bam_files(
      pipeline_results$dedup_output,
      output_dir = config$output_dirs$dedu_meth,
      threads = config$threads
    )
    
    pipeline_results$sorted_bams <- sorted_bams
  } else if ("dedup" %in% steps && !config$need_deduplication) {
    message("\nSkipping deduplication (not required for RRBS)")
    pipeline_results$dedup_output <- pipeline_results$align_output
    
    # Sort BAM files for methylKit
    sorted_bams <- sort_bam_files(
      pipeline_results$align_output,
      output_dir = config$output_dirs$meth_bam,
      threads = config$threads
    )
    
    pipeline_results$sorted_bams <- sorted_bams
  }
  
  # Step 4: Methylation Extraction
  if ("extract" %in% steps) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 4: Methylation Information Extraction")
    message("=" <- rep("=", 60))
    
    # Use sorted BAMs
    bam_files <- if (!is.null(pipeline_results$sorted_bams)) {
      pipeline_results$sorted_bams
    } else if (!is.null(pipeline_results$dedup_output)) {
      pipeline_results$dedup_output
    } else {
      pipeline_results$align_output
    }
    
    genome_path <- file.path(config$ref_dir, toupper(config$genome))
    
    pipeline_results$extract_output <- extract_methylation(
      bam_files,
      genome_dir = genome_path,
      output_dir = config$output_dirs$methinfo,
      config = config
    )
    
    # Generate methylation statistics
    pipeline_results$meth_stats <- calculate_methylation_stats(
      pipeline_results$extract_output
    )
    
    # Print statistics
    print(pipeline_results$meth_stats)
    
    # Generate Bismark summary report
    generate_bismark_report(config$output_dirs$methinfo)
  }
  
  # Step 5: Differential Analysis
  if ("diff" %in% steps) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 5: Differential Methylation Analysis")
    message("=" <- rep("=", 60))
    
    # Check if sample_info is provided
    if (is.null(sample_info)) {
      stop("sample_info is required for differential analysis.")
    }
    
    # Use sorted BAMs
    bam_files <- if (!is.null(pipeline_results$sorted_bams)) {
      pipeline_results$sorted_bams
    } else if (!is.null(pipeline_results$dedup_output)) {
      pipeline_results$dedup_output
    } else {
      pipeline_results$align_output
    }
    
    pipeline_results$diff_output <- run_diff_analysis(
      bam_files,
      sample_info,
      config
    )
  }
  
  # Step 6: Annotation
  if ("annotate" %in% steps) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 6: Annotation and Enrichment")
    message("=" <- rep("=", 60))
    
    if (is.null(pipeline_results$diff_output)) {
      stop("Differential analysis results not found. Please run 'diff' step first.")
    }
    
    # Get DMR data
    dmr_data <- pipeline_results$diff_output$methyl_diff
    
    pipeline_results$annotate_output <- run_annotation(
      dmr_data,
      config
    )
  }
  
  # Step 7: Visualization
  if ("viz" %in% steps) {
    message("\n")
    message("=" <- rep("=", 60))
    message("STEP 7: Visualization")
    message("=" <- rep("=", 60))
    
    if (is.null(pipeline_results$diff_output)) {
      stop("Differential analysis results not found. Please run 'diff' step first.")
    }
    
    # Get DMR data
    dmr_data <- pipeline_results$diff_output$methyl_diff
    
    # Create summary plots
    create_summary_plots(
      dmr_data,
      config,
      output_dir = file.path(config$output_dirs$results, "plots")
    )
  }
  
  # Save pipeline results
  results_file <- file.path(config$output_dirs$results, "pipeline_results.rds")
  saveRDS(pipeline_results, results_file)
  
  message("\n")
  message("=" <- rep("=", 60))
  message("PIPELINE COMPLETED")
  message("=" <- rep("=", 60))
  message("Results saved to:", results_file)
  message("\nOutput directories:")
  for (dir_name in names(config$output_dirs)) {
    message("  ", dir_name, ":", config$output_dirs[[dir_name]])
  }
  
  return(pipeline_results)
}

#' Sort BAM files with samtools
#'
#' @param bam_files Named list of BAM file paths
#' @param output_dir Output directory for sorted BAMs
#' @param threads Number of threads
#' @return List of sorted BAM file paths
#' @export
sort_bam_files <- function(bam_files, output_dir, threads = 4) {
  
  # Check for samtools
  samtools_path <- Sys.which("samtools")
  if (samtools_path == "") {
    stop("samtools not found. Please install: conda install -c bioconda samtools")
  }
  
  message("Sorting BAM files with samtools...")
  
  sorted_bams <- list()
  
  for (sample_name in names(bam_files)) {
    bam_file <- bam_files[[sample_name]]
    
    if (!file.exists(bam_file)) {
      warning("BAM file not found: ", bam_file)
      next
    }
    
    message("  Sorting:", sample_name)
    
    # Output file
    sorted_bam <- file.path(output_dir, paste0(sample_name, "_sorted.bam"))
    
    # Run samtools sort
    cmd <- paste(samtools_path, 
                "sort",
                "-@", threads,
                bam_file,
                "-o", sorted_bam)
    
    result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (result != 0) {
      warning("Sorting failed for sample: ", sample_name)
    } else {
      # Index sorted BAM
      index_cmd <- paste(samtools_path, "index", sorted_bam)
      system(index_cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
      
      sorted_bams[[sample_name]] <- sorted_bam
    }
  }
  
  message("Sorting completed for", length(sorted_bams), "samples.")
  
  return(sorted_bams)
}

#' Generate pipeline summary report
#'
#' @param pipeline_results Results from run_pipeline
#' @param config Configuration list
#' @param output_file Path to save summary report
#' @export
generate_pipeline_summary <- function(pipeline_results, config, 
                                    output_file = "outputs/pipeline_summary.txt") {
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  sink(output_file)
  
  cat("=" <- rep("=", 80), "\n")
  cat("Methylation Analysis Pipeline Summary Report\n")
  cat("=" <- rep("=", 80), "\n\n")
  
  cat("Analysis Type:", config$analysis_type, "\n")
  cat("Genome:", config$genome, "\n")
  cat("Date:", Sys.time(), "\n\n")
  
  # Quality Control
  if (!is.null(pipeline_results$trim_output)) {
    cat("Quality Control:\n")
    cat("  Samples processed:", length(pipeline_results$trim_output), "\n\n")
  }
  
  # Alignment
  if (!is.null(pipeline_results$align_output)) {
    cat("Alignment:\n")
    cat("  Samples aligned:", length(pipeline_results$align_output), "\n\n")
  }
  
  # Deduplication
  if (!is.null(pipeline_results$dedup_output)) {
    cat("Deduplication:\n")
    cat("  Samples deduplicated:", length(pipeline_results$dedup_output), "\n\n")
  }
  
  # Methylation Statistics
  if (!is.null(pipeline_results$meth_stats)) {
    cat("Methylation Statistics:\n")
    print(pipeline_results$meth_stats)
    cat("\n")
  }
  
  # Differential Analysis
  if (!is.null(pipeline_results$diff_output)) {
    cat("Differential Analysis:\n")
    cat("  Hyper-methylated sites:", nrow(pipeline_results$diff_output$hyper_dms), "\n")
    cat("  Hypo-methylated sites:", nrow(pipeline_results$diff_output$hypo_dms), "\n\n")
  }
  
  # Annotation
  if (!is.null(pipeline_results$annotate_output)) {
    cat("Annotation:\n")
    if (!is.null(pipeline_results$annotate_output$go_result)) {
      cat("  GO enriched terms:", nrow(pipeline_results$annotate_output$go_result), "\n")
    }
    if (!is.null(pipeline_results$annotate_output$kegg_result)) {
      cat("  KEGG enriched pathways:", nrow(pipeline_results$annotate_output$kegg_result), "\n")
    }
    cat("\n")
  }
  
  cat("=" <- rep("=", 80), "\n")
  
  sink()
  
  message("Pipeline summary saved to:", output_file)
  
  return(invisible(NULL))
}
