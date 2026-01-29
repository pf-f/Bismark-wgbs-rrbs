#!/usr/bin/env Rscript

#' Install R Dependencies for Methylation Analysis
#'
#' This script installs all required R packages for methylation analysis
#' including Bioconductor packages and CRAN packages.
#'
#' @author Methylation Analysis Toolkit

cat("=" <- rep("=", 70), "\n")
cat("Installing R Dependencies\n")
cat("=" <- rep("=", 70), "\n\n")

# Check for BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(BiocManager))

# Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
bioc_packages <- c(
  "methylKit",
  "genomation", 
  "ChIPseeker",
  "clusterProfiler",
  "AnnotationDbi",
  "GenomicFeatures",
  "GenomicRanges",
  "data.table"
)

for (pkg in bioc_packages) {
  cat("  Installing:", pkg, "...\n")
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# CRAN packages
cat("\nInstalling CRAN packages...\n")
cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "scales",
  "RColorBrewer",
  "viridis",
  "optparse",
  "yaml"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

# Organism databases
cat("\nInstalling organism databases...\n")

# Check for human genome
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat("  Installing: org.Hs.eg.db (human)...\n")
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
} else {
  cat("  Already installed: org.Hs.eg.db\n")
}

# Check for mouse genome
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  cat("  Installing: org.Mm.eg.db (mouse)...\n")
  BiocManager::install("org.Mm.eg.db", ask = FALSE, update = FALSE)
} else {
  cat("  Already installed: org.Mm.eg.db\n")
}

# TxDb packages
cat("\nInstalling TxDb packages (this may take a while)...\n")

# Human hg38
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE)) {
  cat("  Installing: TxDb.Hsapiens.UCSC.hg38.refGene...\n")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene", ask = FALSE, update = FALSE)
} else {
  cat("  Already installed: TxDb.Hsapiens.UCSC.hg38.refGene\n")
}

# Human hg19
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
  cat("  Installing: TxDb.Hsapiens.UCSC.hg19.knownGene...\n")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", ask = FALSE, update = FALSE)
} else {
  cat("  Already installed: TxDb.Hsapiens.UCSC.hg19.knownGene\n")
}

# Mouse mm10
if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
  cat("  Installing: TxDb.Mmusculus.UCSC.mm10.knownGene...\n")
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", ask = FALSE, update = FALSE)
} else {
  cat("  Already installed: TxDb.Mmusculus.UCSC.mm10.knownGene\n")
}

# Optional packages (for advanced visualization)
cat("\nInstalling optional packages (visualization)...\n")

# RIdeogram
if (!requireNamespace("RIdeogram", quietly = TRUE)) {
  cat("  Installing: RIdeogram (for ideogram plots)...\n")
  cat("  Note: This requires installation from GitHub\n")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("jokergoo/RIdeogram")
}

# gggenes
if (!requireNamespace("gggenes", quietly = TRUE)) {
  cat("  Installing: gggenes (for gene visualization)...\n")
  devtools::install_github("wilkox/gggenes")
}

# ggrepel
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  cat("  Installing: ggrepel (for label placement)...\n")
  install.packages("ggrepel", repos = "https://cloud.r-project.org")
}

# Check installation
cat("\n")
cat("=" <- rep("=", 70), "\n")
cat("Installation Summary\n")
cat("=" <- rep("=", 70), "\n\n")

# Check critical packages
critical_packages <- c("methylKit", "genomation", "ChIPseeker", "clusterProfiler")
cat("Critical packages:\n")
for (pkg in critical_packages) {
  status <- if (requireNamespace(pkg, quietly = TRUE)) "OK" else "MISSING"
  cat("  ", pkg, ":", status, "\n")
}

# Check optional packages
optional_packages <- c("RIdeogram", "gggenes", "ggrepel")
cat("\nOptional packages:\n")
for (pkg in optional_packages) {
  status <- if (requireNamespace(pkg, quietly = TRUE)) "OK" else "NOT INSTALLED"
  cat("  ", pkg, ":", status, "\n")
}

cat("\n")
cat("=" <- rep("=", 70), "\n")
cat("Installation completed!\n")
cat("=" <- rep("=", 70), "\n")
cat("\nNote: Some packages may require additional system libraries.\n")
cat("If you encounter errors, please check the package documentation.\n")

# Print R session info
cat("\nR Session Info:\n")
sessionInfo()
