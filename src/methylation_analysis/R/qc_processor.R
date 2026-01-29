#' Quality Control and Preprocessing Module
#'
#' This module handles QC and preprocessing for WGBS/RRBS data
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Run FastQC on raw fastq files
#'
#' @param input_files Character vector of input fastq files
#' @param output_dir Output directory for FastQC results
#' @param threads Number of threads to use
#' @return Path to FastQC output directory
#' @export
run_fastqc <- function(input_files, output_dir = "outputs/fastqc", threads = 4) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if fastqc is available
  fastqc_path <- Sys.which("fastqc")
  if (fastqc_path == "") {
    stop("fastqc not found in PATH. Please install via conda: conda install -c bioconda fastqc")
  }
  
  message("Running FastQC on ", length(input_files), " files...")
  
  # Build command
  cmd <- paste(
    fastqc_path,
    "-t", threads,
    "-o", output_dir,
    paste(input_files, collapse = " ")
  )
  
  # Execute command
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    stop("FastQC failed with exit code: ", result)
  }
  
  message("FastQC completed successfully. Results in: ", output_dir)
  
  return(output_dir)
}

#' Run Trimmomatic for read trimming
#'
#' @param input_files Character vector of paired-end fastq files
#' @param output_dir Output directory for trimmed files
#' @param config Configuration list
#' @return List of output file paths
#' @export
run_trimmomatic <- function(input_files, output_dir, config) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if trimmomatic is available
  trimmomatic_jar <- system2("find", c("/usr", "-name", "trimmomatic.jar", 
                                       "-type", "f", "2>/dev/null"), 
                              stdout = TRUE)
  
  if (length(trimmomatic_jar) == 0) {
    # Try conda path
    conda_prefix <- Sys.getenv("CONDA_PREFIX")
    if (conda_prefix != "") {
      trimmomatic_jar <- list.files(
        file.path(conda_prefix, "share"),
        pattern = "trimmomatic-*.jar",
        full.names = TRUE,
        recursive = TRUE
      )
    }
    
    if (length(trimmomatic_jar) == 0) {
      stop("trimmomatic.jar not found. Please install via conda: conda install -c bioconda trimmomatic")
    }
  }
  
  message("Running Trimmomatic...")
  
  # Java settings
  java_opts <- paste("-Xms20G", "-Xmx20G")
  
  # Trimmomatic parameters
  sliding_window <- config$trim_sliding_window
  min_length <- config$min_read_length
  
  # Process each paired sample
  output_files <- list()
  
  for (r1_file in input_files[grep("_1\\.fastq", input_files)]) {
    # Get R2 file
    sample_name <- gsub("_1\\.fastq.*$", "", basename(r1_file))
    r2_file <- gsub("_1\\.", "_2.", r1_file)
    
    if (!file.exists(r2_file)) {
      warning("R2 file not found for: ", r1_file, ". Skipping.")
      next
    }
    
    # Output files
    r1_output <- file.path(output_dir, paste0(sample_name, "_clip.1.fq.gz"))
    r2_output <- file.path(output_dir, paste0(sample_name, "_clip.2.fq.gz"))
    r1_single <- file.path(output_dir, paste0(sample_name, "_single.R1.fastq.gz"))
    r2_single <- file.path(output_dir, paste0(sample_name, "_single.R2.fastq.gz"))
    
    # Build command
    cmd <- paste(
      "java", java_opts,
      "-jar", trimmomatic_jar[1],
      "PE",
      "-threads", config$threads,
      "-phred33",
      r1_file, r2_file,
      r1_output, r1_single,
      r2_output, r2_single,
      paste0("SLIDINGWINDOW:", sliding_window),
      paste0("MINLEN:", min_length)
    )
    
    message("Processing sample: ", sample_name)
    result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (result != 0) {
      warning("Trimmomatic failed for sample: ", sample_name)
    } else {
      message("Completed: ", sample_name)
      output_files[[sample_name]] <- list(
        r1 = r1_output,
        r2 = r2_output
      )
    }
  }
  
  # Remove single reads (not needed for downstream analysis)
  for (file in list.files(output_dir, pattern = "*single*\\.gz", full.names = TRUE)) {
    file.remove(file)
  }
  
  message("Trimmomatic completed. Processed ", length(output_files), " samples.")
  
  return(output_files)
}

#' Run Trim Galore (recommended for RRBS)
#'
#' @param input_files Character vector of paired-end fastq files
#' @param output_dir Output directory for trimmed files
#' @param config Configuration list
#' @return List of output file paths
#' @export
run_trim_galore <- function(input_files, output_dir, config) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if trim_galore is available
  trim_galore_path <- Sys.which("trim_galore")
  if (trim_galore_path == "") {
    stop("trim_galore not found in PATH. Please install via conda: conda install -c bioconda trim-galore")
  }
  
  message("Running Trim Galore...")
  
  # RRBS-specific parameters
  rrbs_params <- if (config$use_trim_galore && config$trim_galore_rrbs) {
    "--rrbs"
  } else {
    ""
  }
  
  # Process each paired sample
  output_files <- list()
  
  for (r1_file in input_files[grep("_1\\.fastq", input_files)]) {
    # Get R2 file
    r2_file <- gsub("_1\\.", "_2.", r1_file)
    
    if (!file.exists(r2_file)) {
      warning("R2 file not found for: ", r1_file, ". Skipping.")
      next
    }
    
    sample_name <- gsub("_1\\.fastq.*$", "", basename(r1_file))
    
    # Build command
    cmd <- paste(
      trim_galore_path,
      "--paired",
      rrbs_params,
      "--phred33",
      "-o", output_dir,
      r1_file, r2_file
    )
    
    message("Processing sample: ", sample_name)
    result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (result != 0) {
      warning("Trim Galore failed for sample: ", sample_name)
    } else {
      message("Completed: ", sample_name)
      
      # Find output files (trim_galore adds _val_1.fq.gz suffix)
      output_files[[sample_name]] <- list(
        r1 = file.path(output_dir, paste0(sample_name, "_val_1.fq.gz")),
        r2 = file.path(output_dir, paste0(sample_name, "_val_2.fq.gz"))
      )
    }
  }
  
  message("Trim Galore completed. Processed ", length(output_files), " samples.")
  
  return(output_files)
}

#' Main QC function
#'
#' @param input_files Character vector of input fastq files
#' @param output_dir Output directory
#' @param config Configuration list
#' @return List of output file paths
#' @export
run_qc <- function(input_files, output_dir, config) {
  
  message("=" <- rep("=", 60))
  message("Quality Control and Preprocessing")
  message("=" <- rep("=", 60))
  
  # Run FastQC
  fastqc_dir <- file.path(output_dir, "fastqc")
  run_fastqc(input_files, fastqc_dir, config$threads)
  
  # Run trimming
  clip_dir <- config$output_dirs$clip_reads
  
  if (config$use_trim_galore && config$trim_galore_rrbs) {
    output_files <- run_trim_galore(input_files, clip_dir, config)
  } else {
    output_files <- run_trimmomatic(input_files, clip_dir, config)
  }
  
  # Run FastQC again on trimmed files
  if (length(output_files) > 0) {
    trimmed_files <- unlist(output_files)
    fastqc_trim_dir <- file.path(output_dir, "fastqc_trimmed")
    run_fastqc(trimmed_files, fastqc_trim_dir, config$threads)
  }
  
  message("\nQC completed successfully!")
  
  return(output_files)
}

#' Generate MultiQC report
#'
#' @param qc_dir Directory containing QC results
#' @param output_file Output file path for MultiQC report
#' @export
run_multiqc <- function(qc_dir, output_file = "outputs/multiqc_report.html") {
  
  # Check if multiqc is available
  multiqc_path <- Sys.which("multiqc")
  if (multiqc_path == "") {
    message("multiqc not found. Skipping MultiQC report.")
    return(invisible(NULL))
  }
  
  message("Generating MultiQC report...")
  
  cmd <- paste(multiqc_path, qc_dir, "-o", dirname(output_file))
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    warning("MultiQC failed with exit code: ", result)
  } else {
    message("MultiQC report generated: ", output_file)
  }
  
  return(invisible(output_file))
}
