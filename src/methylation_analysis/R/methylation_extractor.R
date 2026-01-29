#' Methylation Information Extraction Module
#'
#' This module handles extraction of methylation information from Bismark BAM files
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Run Bismark methylation extractor
#'
#' @param bam_files Character vector of BAM files
#' @param genome_dir Path to Bismark genome directory
#' @param output_dir Output directory for methylation information
#' @param config Configuration list
#' @return List of output file paths
#' @export
extract_methylation <- function(bam_files, genome_dir, output_dir, config) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if bismark_methylation_extractor is available
  extractor_path <- Sys.which("bismark_methylation_extractor")
  if (extractor_path == "") {
    stop("bismark_methylation_extractor not found. Please ensure bismark is properly installed.")
  }
  
  message("=" <- rep("=", 60))
  message("Extracting Methylation Information")
  message("=" <- rep("=", 60))
  message("Processing ", length(bam_files), " BAM files...")
  
  # Build base command
  cmd_parts <- c(
    extractor_path,
    "--paired-end",
    "--gzip",
    "--parallel", as.character(config$bismark_parallel),
    "--buffer_size", config$bismark_buffer_size,
    "--genome_folder", genome_dir
  )
  
  # Add optional parameters
  if (config$extract_bedgraph) {
    cmd_parts <- c(cmd_parts, "--bedGraph")
  }
  
  if (config$extract_cytosine_report) {
    cmd_parts <- c(cmd_parts, "--cytosine_report", "--report")
  }
  
  if (config$extract_comprehensive) {
    cmd_parts <- c(cmd_parts, "--CX", "--comprehensive")
  }
  
  output_files <- list()
  
  for (bam_file in bam_files) {
    sample_name <- gsub("\\.bam$", "", basename(bam_file))
    
    message("\nExtracting methylation from: ", sample_name)
    
    # Build sample-specific command
    sample_cmd <- c(cmd_parts,
                   bam_file,
                   "-o", output_dir)
    
    cmd <- paste(sample_cmd, collapse = " ")
    
    # Execute command
    start_time <- Sys.time()
    result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    end_time <- Sys.time()
    
    if (result != 0) {
      warning("Methylation extraction failed for sample: ", sample_name)
    } else {
      elapsed <- difftime(end_time, start_time, units = "mins")
      message("  Completed in", round(as.numeric(elapsed), 2), "minutes")
      
      # Collect output files
      output_files[[sample_name]] <- list(
        bam = bam_file,
        cpg = file.path(output_dir, paste0("CpG_context_", sample_name, ".txt.gz")),
        chg = file.path(output_dir, paste0("CHG_context_", sample_name, ".txt.gz")),
        chh = file.path(output_dir, paste0("CHH_context_", sample_name, ".txt.gz"))
      )
      
      # Add optional output files
      if (config$extract_bedgraph) {
        output_files[[sample_name]]$bedgraph <- file.path(
          output_dir, 
          paste0(sample_name, ".bedGraph.gz")
        )
      }
      
      if (config$extract_cytosine_report) {
        output_files[[sample_name]]$cytosine_report <- file.path(
          output_dir, 
          paste0(sample_name, "_CpG_report.txt.gz")
        )
      }
    }
  }
  
  message("\nMethylation extraction completed!")
  message("Successfully processed:", length(output_files), "/", length(bam_files), "samples")
  
  return(output_files)
}

#' Parse Bismark coverage file
#'
#' @param coverage_file Path to Bismark coverage file (uncompressed or gzipped)
#' @param context Methylation context (CpG, CHG, or CHH)
#' @return Data frame with methylation data
#' @export
parse_coverage_file <- function(coverage_file, context = "CpG") {
  
  # Check if file exists
  if (!file.exists(coverage_file)) {
    stop("Coverage file not found: ", coverage_file)
  }
  
  # Read file (handle gzipped files)
  if (grepl("\\.gz$", coverage_file)) {
    con <- gzfile(coverage_file)
  } else {
    con <- file(coverage_file)
  }
  
  data <- tryCatch({
    read.table(con, 
              header = FALSE,
              stringsAsFactors = FALSE,
              sep = "\t",
              comment.char = "#",
              col.names = c("chromosome", "position", "strand",
                           "count_methylated", "count_unmethylated"))
  }, finally = {
    close(con)
  })
  
  if (nrow(data) == 0) {
    warning("No data found in coverage file: ", coverage_file)
    return(data.frame())
  }
  
  # Calculate methylation percentage
  data$total_count <- data$count_methylated + data$count_unmethylated
  data$methylation_percentage <- 
    (data$count_methylated / data$total_count) * 100
  
  # Add context
  data$context <- context
  
  return(data)
}

#' Parse cytosine report
#'
#' @param cytosine_report Path to cytosine report file
#' @return Data frame with cytosine report
#' @export
parse_cytosine_report <- function(cytosine_report) {
  
  # Check if file exists
  if (!file.exists(cytosine_report)) {
    stop("Cytosine report not found: ", cytosine_report)
  }
  
  # Read file
  if (grepl("\\.gz$", cytosine_report)) {
    con <- gzfile(cytosine_report)
  } else {
    con <- file(cytosine_report)
  }
  
  data <- tryCatch({
    read.table(con,
              header = TRUE,
              stringsAsFactors = FALSE,
              sep = "\t")
  }, finally = {
    close(con)
  })
  
  if (nrow(data) == 0) {
    warning("No data found in cytosine report: ", cytosine_report)
    return(data.frame())
  }
  
  return(data)
}

#' Calculate genome-wide methylation statistics
#'
#' @param coverage_files Named list of coverage files
#' @return Data frame with methylation statistics
#' @export
calculate_methylation_stats <- function(coverage_files) {
  
  stats_list <- lapply(names(coverage_files), function(sample_name) {
    files <- coverage_files[[sample_name]]
    
    # Parse each context
    cpg_data <- parse_coverage_file(files$cpg, "CpG")
    chg_data <- parse_coverage_file(files$chg, "CHG")
    chh_data <- parse_coverage_file(files$chh, "CHH")
    
    # Calculate statistics for CpG
    cpg_stats <- if (nrow(cpg_data) > 0) {
      list(
        total_sites = nrow(cpg_data),
        total_reads = sum(cpg_data$total_count),
        methylated_reads = sum(cpg_data$count_methylated),
        avg_methylation = mean(cpg_data$methylation_percentage),
        median_methylation = median(cpg_data$methylation_percentage)
      )
    } else {
      list(
        total_sites = 0,
        total_reads = 0,
        methylated_reads = 0,
        avg_methylation = NA,
        median_methylation = NA
      )
    }
    
    # Calculate statistics for CHG
    chg_stats <- if (nrow(chg_data) > 0) {
      list(
        total_sites = nrow(chg_data),
        avg_methylation = mean(chg_data$methylation_percentage)
      )
    } else {
      list(
        total_sites = 0,
        avg_methylation = NA
      )
    }
    
    # Calculate statistics for CHH
    chh_stats <- if (nrow(chh_data) > 0) {
      list(
        total_sites = nrow(chh_data),
        avg_methylation = mean(chh_data$methylation_percentage)
      )
    } else {
      list(
        total_sites = 0,
        avg_methylation = NA
      )
    }
    
    data.frame(
      Sample = sample_name,
      CpG_Sites = cpg_stats$total_sites,
      CpG_Avg_Meth = round(cpg_stats$avg_methylation, 2),
      CpG_Median_Meth = round(cpg_stats$median_methylation, 2),
      CHG_Sites = chg_stats$total_sites,
      CHG_Avg_Meth = round(chg_stats$avg_methylation, 2),
      CHH_Sites = chh_stats$total_sites,
      CHH_Avg_Meth = round(chh_stats$avg_methylation, 2)
    )
  })
  
  df <- do.call(rbind, stats_list)
  
  return(df)
}

#' Merge coverage files from multiple samples
#'
#' @param coverage_files Named list of coverage file paths
#' @param output_file Output file path for merged data
#' @param context Methylation context (CpG, CHG, or CHH)
#' @export
merge_coverage_files <- function(coverage_files, output_file, context = "CpG") {
  
  # Read and combine all files
  combined_data <- lapply(names(coverage_files), function(sample_name) {
    context_key <- tolower(context)
    if (context_key == "cpg") {
      file_path <- coverage_files[[sample_name]]$cpg
    } else if (context_key == "chg") {
      file_path <- coverage_files[[sample_name]]$chg
    } else if (context_key == "chh") {
      file_path <- coverage_files[[sample_name]]$chh
    } else {
      stop("Invalid context. Must be CpG, CHG, or CHH.")
    }
    
    data <- parse_coverage_file(file_path, context)
    data$sample <- sample_name
    return(data)
  })
  
  merged_data <- do.call(rbind, combined_data)
  
  # Write output
  if (grepl("\\.gz$", output_file)) {
    con <- gzfile(output_file, "w")
  } else {
    con <- file(output_file, "w")
  }
  
  write.table(merged_data, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  
  message("Merged coverage file written to: ", output_file)
  message("Total records:", nrow(merged_data))
  
  return(merged_data)
}

#' Generate Bismark summary report
#'
#' @param output_dir Directory containing Bismark output files
#' @return Path to HTML summary report
#' @export
generate_bismark_report <- function(output_dir) {
  
  # Check if bismark2report is available
  report_path <- Sys.which("bismark2report")
  if (report_path == "") {
    warning("bismark2report not found. Skipping summary report.")
    return(invisible(NULL))
  }
  
  message("Generating Bismark summary report...")
  
  cmd <- paste(report_path, "--dir", output_dir)
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    warning("Bismark report generation failed with exit code: ", result)
  } else {
    # Find generated report
    report_files <- list.files(output_dir, pattern = "*_bismark2report.html$")
    if (length(report_files) > 0) {
      message("Bismark report generated: ", file.path(output_dir, report_files[1]))
      return(file.path(output_dir, report_files[1]))
    }
  }
  
  return(invisible(NULL))
}

#' Filter coverage by depth
#'
#' @param data Data frame from parse_coverage_file
#' @param min_coverage Minimum coverage depth
#' @return Filtered data frame
#' @export
filter_by_coverage <- function(data, min_coverage = 10) {
  
  filtered_data <- data[data$total_count >= min_coverage, ]
  
  message("Filtered", nrow(data), "sites to", nrow(filtered_data), 
          "sites (>= ", min_coverage, "x coverage)")
  
  return(filtered_data)
}
