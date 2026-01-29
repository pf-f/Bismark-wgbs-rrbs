#' Alignment Module using Bismark
#'
#' This module handles genome preparation and alignment with Bismark
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Prepare Bismark genome indices
#'
#' @param genome_dir Path to genome directory containing FASTA files
#' @param aligner_path Path to aligner (bowtie2 or hisat2). If NULL, auto-detect.
#' @param verbose Print verbose output
#' @return Path to Bismark genome directory
#' @export
prepare_bismark_genome <- function(genome_dir, aligner_path = NULL, verbose = TRUE) {
  
  # Check if genome directory exists
  if (!dir.exists(genome_dir)) {
    stop("Genome directory not found: ", genome_dir)
  }
  
  # Check for Bismark genome directory
  bismark_genome_dir <- file.path(genome_dir, "Bisulfite_Genome")
  
  if (dir.exists(bismark_genome_dir)) {
    message("Bismark genome already prepared at: ", bismark_genome_dir)
    message("Skipping genome preparation.")
    return(bismark_genome_dir)
  }
  
  # Check if bismark_genome_preparation is available
  bismark_prep <- Sys.which("bismark_genome_preparation")
  if (bismark_prep == "") {
    stop("bismark_genome_preparation not found. Please install via conda: conda install -c bioconda bismark")
  }
  
  # Check for FASTA files
  fasta_files <- list.files(genome_dir, pattern = "\\.fa$|\\.fasta$", 
                          ignore.case = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No FASTA files found in genome directory: ", genome_dir)
  }
  
  message("Preparing Bismark genome from: ", genome_dir)
  message("Found ", length(fasta_files), " FASTA files")
  
  # Build command
  cmd <- bismark_prep
  if (!is.null(aligner_path)) {
    cmd <- paste(cmd, "--path_to_aligner", aligner_path)
  }
  if (verbose) {
    cmd <- paste(cmd, "--verbose")
  }
  cmd <- paste(cmd, genome_dir)
  
  message("This may take several hours depending on genome size...")
  message("Command: ", cmd)
  
  # Execute command
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    stop("Bismark genome preparation failed with exit code: ", result)
  }
  
  # Verify output
  if (!dir.exists(bismark_genome_dir)) {
    stop("Bismark genome preparation completed but expected directory not found: ", 
         bismark_genome_dir)
  }
  
  message("Bismark genome preparation completed successfully!")
  message("Genome indices at: ", bismark_genome_dir)
  
  return(bismark_genome_dir)
}

#' Run Bismark alignment
#'
#' @param input_files List of paired-end fastq files (output from QC module)
#' @param genome_dir Path to Bismark genome directory
#' @param output_dir Output directory for BAM files
#' @param config Configuration list
#' @return List of output BAM file paths
#' @export
run_bismark <- function(input_files, genome_dir, output_dir, config) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if bismark is available
  bismark_path <- Sys.which("bismark")
  if (bismark_path == "") {
    stop("bismark not found. Please install via conda: conda install -c bioconda bismark")
  }
  
  message("=" <- rep("=", 60))
  message("Running Bismark Alignment")
  message("=" <- rep("=", 60))
  message("Analysis type:", config$analysis_type)
  message("Genome:", config$genome)
  message("Processing ", length(input_files), " samples...")
  
  # Build base command
  cmd_parts <- c(
    bismark_path,
    "--bam",
    "--parallel", as.character(config$bismark_parallel),
    "--fastq",
    "--genome", genome_dir
  )
  
  # RRBS-specific parameters
  if (config$analysis_type == "RRBS") {
    # Bismark automatically handles RRBS data differently
    message("Using RRBS-specific alignment parameters")
  }
  
  output_bams <- list()
  
  # Process each sample
  for (sample_name in names(input_files)) {
    files <- input_files[[sample_name]]
    
    if (is.null(files$r2) || !file.exists(files$r2)) {
      warning("Paired-end files not found for sample: ", sample_name, ". Skipping.")
      next
    }
    
    message("\nAligning sample: ", sample_name)
    message("  R1:", files$r1)
    message("  R2:", files$r2)
    
    # Build sample-specific command
    sample_cmd <- c(cmd_parts,
                   "-1", files$r1,
                   "-2", files$r2,
                   "-o", output_dir,
                   "--prefix", sample_name)
    
    cmd <- paste(sample_cmd, collapse = " ")
    
    # Execute command
    start_time <- Sys.time()
    result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    end_time <- Sys.time()
    
    if (result != 0) {
      warning("Bismark alignment failed for sample: ", sample_name)
    } else {
      elapsed <- difftime(end_time, start_time, units = "mins")
      message("  Completed in", round(as.numeric(elapsed), 2), "minutes")
      
      # Find output BAM file
      bam_pattern <- paste0(sample_name, ".*\\.bam$")
      bam_files <- list.files(output_dir, pattern = bam_pattern, full.names = TRUE)
      
      if (length(bam_files) > 0) {
        output_bams[[sample_name]] <- bam_files[1]
        message("  Output BAM:", bam_files[1])
      }
    }
  }
  
  message("\nBismark alignment completed!")
  message("Successfully aligned:", length(output_bams), "/", length(input_files), "samples")
  
  return(output_bams)
}

#' Check Bismark installation
#'
#' @return Named list with Bismark version and aligner information
#' @export
check_bismark <- function() {
  
  result <- list(
    bismark_installed = FALSE,
    bismark_version = NA,
    aligner = NA,
    aligner_version = NA
  )
  
  # Check Bismark
  bismark_path <- Sys.which("bismark")
  if (bismark_path != "") {
    result$bismark_installed <- TRUE
    
    # Get version
    version_output <- system2(bismark_path, "--version", stdout = TRUE, stderr = TRUE)
    if (length(version_output) > 0) {
      result$bismark_version <- version_output[1]
    }
  }
  
  # Check aligner (prefer bowtie2)
  bowtie2_path <- Sys.which("bowtie2")
  if (bowtie2_path != "") {
    result$aligner <- "bowtie2"
    version_output <- system2(bowtie2_path, "--version", stdout = TRUE)
    if (length(version_output) > 0) {
      result$aligner_version <- version_output[1]
    }
  } else {
    # Check for hisat2
    hisat2_path <- Sys.which("hisat2")
    if (hisat2_path != "") {
      result$aligner <- "hisat2"
      version_output <- system2(hisat2_path, "--version", stdout = TRUE)
      if (length(version_output) > 0) {
        result$aligner_version <- version_output[1]
      }
    }
  }
  
  return(result)
}

#' Print Bismark installation status
#'
#' @export
print_bismark_status <- function() {
  
  status <- check_bismark()
  
  cat("=" <- rep("=", 60), "\n")
  cat("Bismark Installation Status\n")
  cat("=" <- rep("=", 60), "\n\n")
  
  if (status$bismark_installed) {
    cat("Bismark: Installed\n")
    cat("Version:", status$bismark_version, "\n")
  } else {
    cat("Bismark: NOT INSTALLED\n")
    cat("Please install: conda install -c bioconda bismark\n")
  }
  
  if (!is.na(status$aligner)) {
    cat("\nAligner:", status$aligner, "\n")
    cat("Version:", status$aligner_version, "\n")
  } else {
    cat("\nAligner: NOT FOUND\n")
    cat("Please install bowtie2: conda install -c bioconda bowtie2\n")
    cat("Or hisat2: conda install -c bioconda hisat2\n")
  }
  
  cat("\n", "=" <- rep("=", 60), "\n")
  
  invisible(status)
}
