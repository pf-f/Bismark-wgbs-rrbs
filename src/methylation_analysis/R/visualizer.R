#' Visualization Module
#'
#' This module handles visualization of methylation analysis results
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Load required libraries for visualization
#'
#' @export
load_viz_libraries <- function() {
  
  required_packages <- c("ggplot2", "dplyr", "tidyr", "scales",
                        "RColorBrewer", "viridis")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' not installed. ",
           "Please install: install.packages('", pkg, "')")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Optional libraries
  optional_packages <- list(
    "RIdeogram" = list(
      install = function() {
        if (!requireNamespace("AnnoProbe", quietly = TRUE)) {
          install.packages("AnnoProbe")
        }
        install_github("jokergoo/ComplexHeatmap")
      },
      url = "https://github.com/jokergoo/RIdeogram"
    ),
    "gggenes" = list(
      install = function() {
        install.packages("gggenes")
      },
      url = "https://github.com/wilkox/gggenes"
    ),
    "ggrepel" = list(
      install = function() {
        install.packages("ggrepel")
      },
      url = "https://github.com/slowkow/ggrepel"
    )
  )
  
  # Try to load optional libraries
  for (pkg in names(optional_packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Optional package '", pkg, "' not installed.")
      message("  Install for advanced visualization features:")
      message("  ", optional_packages[[pkg]]$url)
    }
  }
  
  message("Loaded visualization libraries.")
  
  invisible(TRUE)
}

#' Plot methylation statistics
#'
#' @param methyl_stats Data frame with methylation statistics
#' @param output_file Path to save plot
#' @export
plot_methylation_statistics <- function(methyl_stats, 
                                       output_file = "outputs/plots/methylation_stats.png") {
  
  load_viz_libraries()
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Convert to long format
  stats_long <- methyl_stats %>%
    select(Sample, CpG_Avg_Meth, CHG_Avg_Meth, CHH_Avg_Meth) %>%
    pivot_longer(cols = -Sample, 
                names_to = "Context", 
                values_to = "Methylation")
  
  # Create plot
  p <- ggplot(stats_long, aes(x = Sample, y = Methylation, fill = Context)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Average Methylation by Context",
         x = "Sample", 
         y = "Methylation Percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  
  message("Saved methylation statistics plot:", output_file)
  
  return(invisible(p))
}

#' Plot DMR chromosome distribution
#'
#' @param dmr_data Data frame with DMR information
#' @param output_file Path to save plot
#' @export
plot_dmr_chromosome_distribution <- function(dmr_data, 
                                            output_file = "outputs/plots/dmr_chr_dist.png") {
  
  load_viz_libraries()
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Count DMRs per chromosome
  chr_counts <- dmr_data %>%
    group_by(chr) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  # Ensure correct chromosome order
  chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
  chr_counts$chr <- factor(chr_counts$chr, levels = chr_order)
  chr_counts <- chr_counts[!is.na(chr_counts$chr), ]
  
  # Create plot
  p <- ggplot(chr_counts, aes(x = chr, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "DMR Distribution Across Chromosomes",
         x = "Chromosome", 
         y = "Number of DMRs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  
  message("Saved chromosome distribution plot:", output_file)
  
  return(invisible(p))
}

#' Plot DMRs on chromosomes using ggplot2
#'
#' @param dmr_data Data frame with DMR information
#' @param genome Genome version (e.g., "hg38", "mm10")
#' @param output_file Path to save plot
#' @export
plot_dmr_genome <- function(dmr_data, genome = "hg38",
                             output_file = "outputs/plots/dmr_genome.png") {
  
  load_viz_libraries()
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Get chromosome lengths
  chrom_lengths <- get_chromosome_lengths(genome)
  
  # Ensure chromosome order
  chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")
  
  # Filter DMRs to chromosomes in reference
  dmr_data <- dmr_data[dmr_data$chr %in% chrom_order, ]
  dmr_data$chr <- factor(dmr_data$chr, levels = chrom_order)
  
  # Create y-axis positions
  y_positions <- setNames(rev(1:length(chrom_order)), chrom_order)
  
  # Add y position to DMR data
  dmr_data$y_pos <- y_positions[as.character(dmr_data$chr)]
  
  # Create plot
  p <- ggplot() +
    # Draw chromosomes
    geom_rect(data = chrom_lengths,
              aes(xmin = 0, xmax = length,
                   ymin = y_positions[as.character(chr)] - 0.4,
                   ymax = y_positions[as.character(chr)] + 0.4),
              fill = "lightgray", color = "black") +
    # Draw DMRs
    geom_point(data = dmr_data,
              aes(x = (start + end) / 2, 
                   y = y_pos,
                   color = meth.diff,
                   size = abs(meth.diff)),
              alpha = 0.7) +
    # Set color gradient
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, 
                          name = "Methylation\nDifference") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_y_continuous(breaks = 1:length(y_positions),
                       labels = names(y_positions)) +
    scale_x_continuous(labels = scales::comma_format(scale = 1e-6, suffix = "M")) +
    labs(title = paste("DMR Distribution on", toupper(genome), "Genome"),
         subtitle = "Points represent DMR regions, colored by methylation difference",
         x = "Genomic Position (Mb)",
         y = "Chromosome") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right")
  
  ggsave(output_file, p, width = 14, height = 8, dpi = 300)
  
  message("Saved genome plot:", output_file)
  
  return(invisible(p))
}

#' Plot DMR ideogram (requires RIdeogram)
#'
#' @param dmr_data Data frame with DMR information
#' @param output_file Path to save plot
#' @export
plot_dmr_ideogram <- function(dmr_data, 
                               output_file = "outputs/plots/dmr_ideogram.svg") {
  
  # Check for RIdeogram
  if (!requireNamespace("RIdeogram", quietly = TRUE)) {
    message("RIdeogram not installed. Skipping ideogram plot.")
    message("  Install: install_github('jokergoo/RIdeogram')")
    return(invisible(NULL))
  }
  
  library(RIdeogram)
  library(dplyr)
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Prepare data for RIdeogram
  # Convert meth.diff to shape
  plot_data <- dmr_data %>%
    select(chr, start, end, meth.diff) %>%
    mutate(
      shape = ifelse(meth.diff > 0, "triangle", "circle"),
      color = meth.diff
    )
  
  # Generate color gradient
  color_grad <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Normalize meth.diff to 1-100
  meth_norm <- scales::rescale(plot_data$color, to = c(1, 100))
  plot_data$color <- color_grad[round(meth_norm)]
  
  # Load karyotype data
  data(human_karyotype, package = "RIdeogram")
  
  # Generate ideogram
  message("Generating RIdeogram plot...")
  
  ideogram(
    karyotype = human_karyotype,
    overlaid = NULL,
    label = plot_data,
    label_type = "marker",
    output = output_file
  )
  
  # Convert to PNG
  if (requireNamespace("svg2png", quietly = TRUE)) {
    png_file <- sub("\\.svg$", ".png", output_file)
    svg2png::svg2png(output_file, png_file)
    message("Saved ideogram (PNG):", png_file)
  }
  
  message("Saved ideogram (SVG):", output_file)
  
  return(invisible(NULL))
}

#' Plot volcano plot for DMRs
#'
#' @param dmr_data Data frame with DMR information
#' @param output_file Path to save plot
#' @export
plot_dmr_volcano <- function(dmr_data, 
                             output_file = "outputs/plots/dmr_volcano.png") {
  
  load_viz_libraries()
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Add significance column
  dmr_data$significant <- dmr_data$qvalue < 0.05
  
  # Create volcano plot
  p <- ggplot(dmr_data, aes(x = meth.diff, y = -log10(qvalue))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                      labels = c("FALSE" = "Not Significant", "TRUE" = "Significant")) +
    geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = "Volcano Plot of DMRs",
         subtitle = "Red points are significant (q-value < 0.05, |diff| > 25%)",
         x = "Methylation Difference (%)",
         y = "-log10(Q-value)",
         color = "Significance") +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(output_file, p, width = 10, height = 8, dpi = 300)
  
  message("Saved volcano plot:", output_file)
  
  return(invisible(p))
}

#' Plot heatmap of DMR methylation
#'
#' @param dmr_data Data frame with DMR information
#' @param sample_info Data frame with sample information
#' @param output_file Path to save plot
#' @export
plot_dmr_heatmap <- function(dmr_data, sample_info,
                             output_file = "outputs/plots/dmr_heatmap.png") {
  
  load_viz_libraries()
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Create heatmap matrix from methyl.diff
  # This is a simplified version
  # Full implementation would require per-sample methylation data
  
  message("Creating DMR heatmap...")
  
  # For now, create a simple summary plot
  p <- ggplot(dmr_data, aes(x = chr, fill = meth.diff)) +
    geom_histogram(bins = 50) +
    facet_wrap(~ chr, scales = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(title = "DMR Methylation Distribution by Chromosome",
         x = "Chromosome",
         y = "Count",
         fill = "Methylation Difference") +
    theme_minimal()
  
  ggsave(output_file, p, width = 14, height = 10, dpi = 300)
  
  message("Saved heatmap:", output_file)
  
  return(invisible(p))
}

#' Get chromosome lengths for genome
#'
#' @param genome Genome version
#' @return Data frame with chromosome lengths
#' @export
get_chromosome_lengths <- function(genome = "hg38") {
  
  lengths <- switch(genome,
                   "hg38" = data.frame(
                     chr = c(paste0("chr", 1:22), "chrX", "chrY"),
                     length = c(248956422, 242193529, 198295559, 190214555, 
                                181538259, 170805979, 159345973, 145138636,
                                138394717, 133797422, 135086622, 133275309,
                                114364328, 107043718, 101991189, 90338345,
                                83257441, 80373285, 58617616, 64444167,
                                46709983, 50818468, 156040895, 57227415)
                   ),
                   "mm10" = data.frame(
                     chr = c(paste0("chr", 1:19), "chrX", "chrY"),
                     length = c(195471971, 182113224, 160039680, 156508116,
                                151834684, 149736546, 145441459, 129401213,
                                124595110, 130694993, 122082543, 120129022,
                                120421639, 124902244, 104043685, 98207768,
                                94987271, 90702639, 61431566, 171031299,
                                91744698)
                   ))
  
  if (is.null(lengths)) {
    stop("Chromosome lengths not available for genome: ", genome)
  }
  
  return(lengths)
}

#' Create summary report with all plots
#'
#' @param dmr_data Data frame with DMR information
#' @param config Configuration list
#' @param output_dir Directory to save all plots
#' @export
create_summary_plots <- function(dmr_data, config, output_dir = "outputs/plots") {
  
  load_viz_libraries()
  
  message("\n")
  message("=" <- rep("=", 60))
  message("Creating Visualization Plots")
  message("=" <- rep("=", 60))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create plots
  plot_dmr_chromosome_distribution(
    dmr_data,
    output_file = file.path(output_dir, "dmr_chr_distribution.png")
  )
  
  plot_dmr_genome(
    dmr_data,
    genome = config$genome,
    output_file = file.path(output_dir, "dmr_genome.png")
  )
  
  plot_dmr_volcano(
    dmr_data,
    output_file = file.path(output_dir, "dmr_volcano.png")
  )
  
  # Try ideogram (optional)
  plot_dmr_ideogram(
    dmr_data,
    output_file = file.path(output_dir, "dmr_ideogram.svg")
  )
  
  message("\nAll plots saved to:", output_dir)
  
  return(invisible(TRUE))
}
