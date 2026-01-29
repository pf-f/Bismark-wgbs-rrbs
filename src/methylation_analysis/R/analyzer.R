#' Differential Methylation Analysis Module
#'
#' This module handles DMR and DMS analysis using methylKit
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Load required libraries
#'
#' @export
load_methylkit <- function() {
  
  required_packages <- c("methylKit", "GenomicRanges", "data.table")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' not installed. Please install: ",
           "BiocManager::install('", pkg, "')")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  message("Loaded methylKit and dependencies.")
  
  invisible(TRUE)
}

#' Process Bismark alignment files with methylKit
#'
#' @param bam_files Named list of BAM file paths
#' @param sample_info Data frame with sample information
#' @param read_context Methylation context (CpG, CHG, or CHH)
#' @param save_folder Folder to save processed data
#' @return methylRawList object
#' @export
process_bismark_for_methylkit <- function(bam_files, sample_info, 
                                          read_context = "CpG", 
                                          save_folder = "outputs/methylkit") {
  
  load_methylkit()
  
  # Create save folder
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  
  # Prepare file list and sample IDs
  file_list <- unlist(bam_files)
  sample_ids <- names(bam_files)
  
  # Match treatment codes from sample_info
  treatment <- sample_info$treatment_code[match(sample_ids, sample_info$sample_alias)]
  
  # Remove any NA values
  valid_idx <- !is.na(treatment)
  file_list <- file_list[valid_idx]
  sample_ids <- sample_ids[valid_idx]
  treatment <- treatment[valid_idx]
  
  message("=" <- rep("=", 60))
  message("Processing Bismark files with methylKit")
  message("=" <- rep("=", 60))
  message("Samples:", length(file_list))
  message("Context:", read_context)
  
  # Process Bismark alignments
  my.methRaw <- tryCatch({
    processBismarkAln(
      location = file_list,
      sample.id = sample_ids,
      treatment = treatment,
      assembly = "hg38",
      read.context = read_context,
      mincov = 0,
      save.folder = save_folder,
      pipeline = "Bismark"
    )
  }, error = function(e) {
    stop("Error processing Bismark files: ", e$message)
  })
  
  message("Successfully processed", length(my.methRaw), "samples")
  
  return(my.methRaw)
}

#' Calculate methylation statistics
#'
#' @param methyl_obj methylRaw or methylBase object
#' @param plot Boolean, whether to generate plots
#' @param save_dir Directory to save plots
#' @return Data frame with statistics
#' @export
get_methylation_stats <- function(methyl_obj, plot = TRUE, save_dir = "outputs/plots") {
  
  load_methylkit()
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # If methylRawList, process each sample
  if (inherits(methyl_obj, "methylRawList")) {
    stats_list <- lapply(seq_along(methyl_obj), function(i) {
      obj <- methyl_obj[[i]]
      
      # Get statistics
      stats <- methylKit::getMethylationStats(obj, plot = FALSE, both.strands = FALSE)
      
      # Generate plot if requested
      if (plot) {
        plot_file <- file.path(save_dir, paste0("methylation_stats_", 
                                               methyl_obj[[i]]@sample.id, ".png"))
        png(plot_file, width = 800, height = 600)
        methylKit::getMethylationStats(obj, plot = TRUE, both.strands = FALSE,
                                      main = paste("Methylation Stats:", 
                                                  methyl_obj[[i]]@sample.id))
        dev.off()
        message("Saved plot:", plot_file)
      }
      
      data.frame(
        Sample = methyl_obj[[i]]@sample.id,
        Sites_Analyzed = stats$totalCs,
        Avg_Coverage = stats$avgCoverge,
        Median_Coverage = stats$medianCoverage
      )
    })
    
    df <- do.call(rbind, stats_list)
  } else {
    df <- methylKit::getMethylationStats(methyl_obj, plot = FALSE, both.strands = FALSE)
  }
  
  return(df)
}

#' Filter by coverage
#'
#' @param methyl_obj methylRawList object
#' @param lo.count Minimum coverage count
#' @param lo.perc Minimum coverage percentile
#' @param hi.count Maximum coverage count
#' @param hi.perc Maximum coverage percentile
#' @return Filtered methylRawList object
#' @export
filter_methylation_coverage <- function(methyl_obj, lo.count = 10, lo.perc = NULL,
                                        hi.count = NULL, hi.perc = 99.9) {
  
  load_methylkit()
  
  message("Filtering by coverage...")
  message("  Minimum coverage:", ifelse(is.null(lo.count), "NULL", lo.count))
  message("  Maximum percentile:", hi.perc, "%")
  
  filtered_obj <- filterByCoverage(
    methyl_obj,
    lo.count = lo.count,
    lo.perc = lo.perc,
    hi.count = hi.count,
    hi.perc = hi.perc
  )
  
  message("Filtering completed.")
  
  return(filtered_obj)
}

#' Unite methylRaw objects
#'
#' @param methyl_obj methylRawList object
#' @param destrand Boolean, whether to destrand (only for CpG)
#' @param min.per.group Minimum samples per group required
#' @return methylBase object
#' @export
unite_methyl_data <- function(methyl_obj, destrand = FALSE, min.per.group = 1) {
  
  load_methylkit()
  
  message("Uniting methylRaw objects...")
  
  united_obj <- tryCatch({
    unite(methyl_obj, 
           destrand = destrand,
           min.per.group = min.per.group)
  }, error = function(e) {
    stop("Error uniting methylRaw objects: ", e$message)
  })
  
  message("Successfully united data.")
  message("  Sites:", nrow(united_obj))
  message("  Samples:", ncol(united_obj) - 3)
  
  return(united_obj)
}

#' Calculate correlation between samples
#'
#' @param methyl_obj methylBase object
#' @param plot Boolean, whether to generate plot
#' @param save_dir Directory to save plot
#' @return Correlation matrix
#' @export
calculate_correlation <- function(methyl_obj, plot = TRUE, save_dir = "outputs/plots") {
  
  load_methylkit()
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  message("Calculating sample correlation...")
  
  # Calculate correlation
  cor_matrix <- getCorrelation(methyl_obj, plot = FALSE)
  
  if (plot) {
    plot_file <- file.path(save_dir, "correlation_heatmap.png")
    png(plot_file, width = 800, height = 600)
    getCorrelation(methyl_obj, plot = TRUE)
    dev.off()
    message("Saved correlation plot:", plot_file)
  }
  
  return(cor_matrix)
}

#' Cluster samples
#'
#' @param methyl_obj methylBase object
#' @param dist Distance method
#' @param method Clustering method
#' @param plot Boolean, whether to generate plot
#' @param save_dir Directory to save plot
#' @return Dendrogram object
#' @export
cluster_samples <- function(methyl_obj, dist = "correlation", method = "ward.D",
                            plot = TRUE, save_dir = "outputs/plots") {
  
  load_methylkit()
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  message("Clustering samples...")
  message("  Distance method:", dist)
  message("  Clustering method:", method)
  
  # Perform clustering
  hc <- clusterSamples(methyl_obj, dist = dist, method = method, plot = FALSE)
  
  if (plot) {
    plot_file <- file.path(save_dir, "sample_clustering.png")
    png(plot_file, width = 1000, height = 600)
    clusterSamples(methyl_obj, dist = dist, method = method, plot = TRUE)
    dev.off()
    message("Saved clustering plot:", plot_file)
  }
  
  return(hc)
}

#' Perform PCA analysis
#'
#' @param methyl_obj methylBase object
#' @param screeplot Boolean, whether to generate screeplot
#' @param plot Boolean, whether to generate PCA plot
#' @param save_dir Directory to save plots
#' @export
perform_pca <- function(methyl_obj, screeplot = TRUE, plot = TRUE, 
                         save_dir = "outputs/plots") {
  
  load_methylkit()
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  message("Performing PCA analysis...")
  
  if (screeplot) {
    plot_file <- file.path(save_dir, "pca_screeplot.png")
    png(plot_file, width = 600, height = 600)
    PCASamples(methyl_obj, screeplot = TRUE)
    dev.off()
    message("Saved PCA screeplot:", plot_file)
  }
  
  if (plot) {
    plot_file <- file.path(save_dir, "pca_plot.png")
    png(plot_file, width = 800, height = 600)
    PCASamples(methyl_obj, screeplot = FALSE)
    dev.off()
    message("Saved PCA plot:", plot_file)
  }
}

#' Calculate differential methylation
#'
#' @param methyl_obj methylBase object
#' @param difference Minimum methylation difference threshold (%)
#' @param qvalue Q-value threshold
#' @param mc.cores Number of cores for parallel processing
#' @return methylDiff object
#' @export
calculate_diff_methylation <- function(methyl_obj, difference = 25, qvalue = 0.01,
                                       mc.cores = 1) {
  
  load_methylkit()
  
  message("=" <- rep("=", 60))
  message("Calculating Differential Methylation")
  message("=" <- rep("=", 60))
  message("  Difference threshold:", difference, "%")
  message("  Q-value threshold:", qvalue)
  message("  Cores:", mc.cores)
  
  # Calculate differential methylation
  myDiff <- tryCatch({
    calculateDiffMeth(methyl_obj, 
                      overdispersion = "MN",
                      test = "Chisq",
                      mc.cores = mc.cores)
  }, error = function(e) {
    stop("Error calculating differential methylation: ", e$message)
  })
  
  message("Differential methylation calculation completed.")
  message("  Total sites tested:", nrow(myDiff))
  
  return(myDiff)
}

#' Get differentially methylated sites/regions
#'
#' @param methyl_diff methylDiff object
#' @param difference Minimum methylation difference (%)
#' @param qvalue Q-value threshold
#' @param type Type of difference: "hyper", "hypo", or "all"
#' @return Data frame with DMS/DMR
#' @export
get_differential_methylated <- function(methyl_diff, difference = 25, qvalue = 0.01,
                                       type = "all") {
  
  load_methylkit()
  
  if (type == "hyper") {
    result <- getMethylDiff(methyl_diff, difference = difference, 
                             qvalue = qvalue, type = "hyper")
  } else if (type == "hypo") {
    result <- getMethylDiff(methyl_diff, difference = difference, 
                             qvalue = qvalue, type = "hypo")
  } else {
    result <- getMethylDiff(methyl_diff, difference = difference, 
                             qvalue = qvalue, type = "all")
  }
  
  message("Found", nrow(result), "differentially methylated sites/regions (", type, ")")
  
  return(result)
}

#' Perform tile-based DMR analysis
#'
#' @param methyl_obj methylBase object
#' @param win.size Window size (bp)
#' @param step.size Step size (bp)
#' @param cov.bases Minimum coverage for bases in window
#' @return methylBase object with tiled data
#' @export
perform_tiling_analysis <- function(methyl_obj, win.size = 1000, step.size = 1000,
                                   cov.bases = 10) {
  
  load_methylkit()
  
  message("Performing tiling analysis...")
  message("  Window size:", win.size, "bp")
  message("  Step size:", step.size, "bp")
  message("  Min coverage:", cov.bases)
  
  tiles <- tileMethylCounts(methyl_obj,
                            win.size = win.size,
                            step.size = step.size,
                            cov.bases = cov.bases)
  
  message("Tiling completed.")
  message("  Tiles generated:", sum(sapply(tiles, nrow)))
  
  return(tiles)
}

#' Run complete differential methylation analysis
#'
#' @param bam_files Named list of BAM file paths
#' @param sample_info Data frame with sample information
#' @param config Configuration list
#' @return List with analysis results
#' @export
run_diff_analysis <- function(bam_files, sample_info, config) {
  
  message("\n")
  message("=" <- rep("=", 60))
  message("Differential Methylation Analysis")
  message("=" <- rep("=", 60))
  
  # Process Bismark files
  save_folder <- file.path(config$output_dirs$results, "methylkit")
  methyl_raw <- process_bismark_for_methylkit(
    bam_files, sample_info, 
    read_context = config$read_context,
    save_folder = save_folder
  )
  
  # Filter by coverage
  methyl_filtered <- filter_methylation_coverage(
    methyl_raw,
    lo.count = config$min_coverage,
    hi.perc = config$max_coverage_percentile
  )
  
  # Unite data
  methyl_base <- unite_methyl_data(methyl_filtered, destrand = FALSE)
  
  # Calculate correlation and clustering
  plots_dir <- file.path(config$output_dirs$results, "plots")
  calculate_correlation(methyl_base, plot = TRUE, save_dir = plots_dir)
  cluster_samples(methyl_base, plot = TRUE, save_dir = plots_dir)
  perform_pca(methyl_base, plot = TRUE, save_dir = plots_dir)
  
  # Calculate differential methylation
  methyl_diff <- calculate_diff_methylation(
    methyl_base,
    difference = config$diff_meth_threshold,
    qvalue = config$diff_qvalue_threshold,
    mc.cores = config$threads
  )
  
  # Get hypermethylated and hypomethylated sites
  hyper_dms <- get_differential_methylated(
    methyl_diff,
    difference = config$diff_meth_threshold,
    qvalue = config$diff_qvalue_threshold,
    type = "hyper"
  )
  
  hypo_dms <- get_differential_methylated(
    methyl_diff,
    difference = config$diff_meth_threshold,
    qvalue = config$diff_qvalue_threshold,
    type = "hypo"
  )
  
  # Return results
  results <- list(
    methyl_raw = methyl_raw,
    methyl_base = methyl_base,
    methyl_diff = methyl_diff,
    hyper_dms = hyper_dms,
    hypo_dms = hypo_dms
  )
  
  # Save results
  results_file <- file.path(save_folder, "diff_analysis_results.rds")
  saveRDS(results, results_file)
  message("\nResults saved to:", results_file)
  
  return(results)
}
