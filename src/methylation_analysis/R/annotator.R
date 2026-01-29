#' Annotation and Enrichment Module
#'
#' This module handles DMR/DMS annotation and functional enrichment
#'
#' @author Methylation Analysis Toolkit
#' @export

#' Load required libraries for annotation
#'
#' @export
load_annotation_libraries <- function() {
  
  required_packages <- c("genomation", "ChIPseeker", "org.Hs.eg.db", 
                        "clusterProfiler", "AnnotationDbi", "GenomicFeatures")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' not installed. Please install: ",
           "BiocManager::install('", pkg, "')")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  message("Loaded annotation libraries.")
  
  invisible(TRUE)
}

#' Read transcript features from BED file
#'
#' @param bed_file Path to BED file with transcript annotations
#' @return List with genomic features
#' @export
read_transcript_features <- function(bed_file) {
  
  load_annotation_libraries()
  
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  
  message("Reading transcript features from:", bed_file)
  
  gene.obj <- readTranscriptFeatures(bed_file)
  
  message("Loaded", length(gene.obj), "features.")
  message("  Features available:", paste(names(gene.obj), collapse = ", "))
  
  return(gene.obj)
}

#' Annotate DMRs/DMSs with gene parts
#'
#' @param dmr_data Data frame with DMR/DMS information
#' @param gene.obj Transcript features object from genomation
#' @return Annotated data frame
#' @export
annotate_with_gene_parts <- function(dmr_data, gene.obj) {
  
  load_annotation_libraries()
  
  # Convert to GRanges
  dmr_gr <- data.frame_to_granges(dmr_data)
  
  message("Annotating with gene parts...")
  
  # Annotate
  annotated <- annotateWithGeneParts(dmr_gr, gene.obj)
  
  # Convert back to data frame
  annotated_df <- as.data.frame(annotated)
  
  message("Annotation completed.")
  
  return(annotated_df)
}

#' Annotate DMRs/DMSs with CpG islands and shores
#'
#' @param dmr_data Data frame with DMR/DMS information
#' @param cpgi_gr GRanges of CpG islands
#' @param shores_gr GRanges of CpG shores
#' @param feature_name Name of the feature
#' @param flank_name Name of the flank
#' @return Annotated GRanges object
#' @export
annotate_with_cpg_features <- function(dmr_data, cpgi_gr, shores_gr,
                                        feature_name = "CpGi", 
                                        flank_name = "shores") {
  
  load_annotation_libraries()
  
  # Convert to GRanges
  dmr_gr <- data.frame_to_granges(dmr_data)
  
  message("Annotating with CpG features...")
  
  # Annotate
  annotated <- annotateWithFeatureFlank(dmr_gr, cpgi_gr, shores_gr,
                                      feature.name = feature_name,
                                      flank.name = flank_name)
  
  message("Annotation completed.")
  
  return(annotated)
}

#' Annotate DMRs/DMSs using ChIPseeker
#'
#' @param dmr_data Data frame with DMR/DMS information
#' @param txdb TxDb object for annotation
#' @param annoDb AnnotationDb for gene symbols
#' @return Annotated data frame
#' @export
annotate_with_chipseeker <- function(dmr_data, txdb, annoDb = NULL) {
  
  load_annotation_libraries()
  
  # Convert to GRanges
  dmr_gr <- data.frame_to_granges(dmr_data)
  
  message("Annotating with ChIPseeker...")
  
  # Annotate
  annotated <- tryCatch({
    annotatePeak(dmr_gr, TxDb = txdb, annoDb = annoDb, verbose = FALSE)
  }, error = function(e) {
    warning("ChIPseeker annotation failed: ", e$message)
    # Return basic annotation without annoDb
    annotatePeak(dmr_gr, TxDb = txdb, verbose = FALSE)
  })
  
  annotated_df <- as.data.frame(annotated)
  
  message("Annotation completed.")
  message("  Annotated regions:", nrow(annotated_df))
  
  return(annotated_df)
}

#' Get genes associated with TSS
#'
#' @param annotated_obj Annotated object from annotateWithGeneParts
#' @param distance Maximum distance to TSS
#' @return Data frame with TSS associations
#' @export
get_tss_association <- function(annotated_obj, distance = 5000) {
  
  load_annotation_libraries()
  
  tss_assoc <- getAssociationWithTSS(annotated_obj, 
                                     d = distance,
                                      getNearest = TRUE)
  
  return(tss_assoc)
}

#' Plot target annotation statistics
#'
#' @param annotated_obj Annotated object
#' @param output_file Path to save plot
#' @export
plot_annotation_stats <- function(annotated_obj, output_file = "outputs/plots/annotation_stats.png") {
  
  load_annotation_libraries()
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  message("Generating annotation statistics plot...")
  
  png(output_file, width = 800, height = 600)
  plotTargetAnnotation(annotated_obj, precedence = TRUE,
                       main = "Differential Methylation Annotation")
  dev.off()
  
  message("Saved plot:", output_file)
  
  return(invisible(NULL))
}

#' Run GO enrichment analysis
#'
#' @param genes Character vector of gene IDs (ENTREZID)
#' @param organism Organism code (e.g., "hsa" for human)
#' @param ontology GO ontology: "BP", "MF", "CC", or "ALL"
#' @param pvalue.cutoff P-value cutoff
#' @param qvalue.cutoff Adjusted p-value cutoff
#' @return enrichResult object
#' @export
run_go_enrichment <- function(genes, organism = "hsa", ontology = "ALL",
                               pvalue.cutoff = 0.05, qvalue.cutoff = 0.2) {
  
  load_annotation_libraries()
  
  # Get appropriate organism database
  org_db <- switch(organism,
                   "hsa" = org.Hs.eg.db,
                   "mmu" = org.Mm.eg.db,
                   "rno" = org.Rn.eg.db,
                   stop("Unsupported organism: ", organism))
  
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Organism database not installed: ", org_db, 
         "\nPlease install: BiocManager::install('", org_db, "')")
  }
  
  message("Running GO enrichment analysis...")
  message("  Genes:", length(genes))
  message("  Ontology:", ontology)
  
  # Run enrichment
  go_result <- tryCatch({
    enrichGO(gene = genes,
              OrgDb = org_db,
              keyType = "ENTREZID",
              ont = ontology,
              pvalueCutoff = pvalue.cutoff,
              qvalueCutoff = qvalue.cutoff,
              minGSSize = 2,
              maxGSSize = 500,
              readable = TRUE)
  }, error = function(e) {
    warning("GO enrichment failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(go_result)) {
    message("GO enrichment completed.")
    message("  Enriched terms:", nrow(go_result))
  }
  
  return(go_result)
}

#' Run KEGG pathway enrichment
#'
#' @param genes Character vector of gene IDs (ENTREZID)
#' @param organism Organism code (e.g., "hsa" for human)
#' @param pvalue.cutoff P-value cutoff
#' @param qvalue.cutoff Adjusted p-value cutoff
#' @return enrichResult object
#' @export
run_kegg_enrichment <- function(genes, organism = "hsa",
                                 pvalue.cutoff = 0.05, qvalue.cutoff = 0.2) {
  
  load_annotation_libraries()
  
  message("Running KEGG enrichment analysis...")
  message("  Genes:", length(genes))
  message("  Organism:", organism)
  
  # Run enrichment
  kegg_result <- tryCatch({
    enrichKEGG(gene = genes,
                organism = organism,
                pvalueCutoff = pvalue.cutoff,
                qvalueCutoff = qvalue.cutoff,
                minGSSize = 2,
                maxGSSize = 500)
  }, error = function(e) {
    warning("KEGG enrichment failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(kegg_result)) {
    message("KEGG enrichment completed.")
    message("  Enriched pathways:", nrow(kegg_result))
  }
  
  return(kegg_result)
}

#' Plot enrichment results
#'
#' @param enrich_result enrichResult object
#' @param plot_type Type of plot: "dot", "bar", "cnet", or "emap"
#' @param output_file Path to save plot
#' @export
plot_enrichment <- function(enrich_result, plot_type = "dot",
                            output_file = "outputs/plots/enrichment.png") {
  
  load_annotation_libraries()
  
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    warning("No enrichment results to plot.")
    return(invisible(NULL))
  }
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  message("Generating enrichment plot (", plot_type, ")...")
  
  png(output_file, width = 1000, height = 600)
  
  result <- switch(plot_type,
                  "dot" = dotplot(enrich_result, showCategory = 20),
                  "bar" = barplot(enrich_result, showCategory = 20),
                  "cnet" = cnetplot(enrich_result, showCategory = 15),
                  "emap" = emapplot(enrich_result, showCategory = 20),
                  stop("Invalid plot_type: ", plot_type))
  
  dev.off()
  
  message("Saved plot:", output_file)
  
  return(invisible(result))
}

#' Convert data frame to GRanges
#'
#' @param df Data frame with columns: chr, start, end (or start/position)
#' @return GRanges object
#' @export
data_frame_to_granges <- function(df) {
  
  load_annotation_libraries()
  
  # Check required columns
  if (!"chr" %in% names(df) && !"chromosome" %in% names(df)) {
    stop("Data frame must have 'chr' or 'chromosome' column")
  }
  
  if (!"start" %in% names(df) && !"position" %in% names(df)) {
    stop("Data frame must have 'start' or 'position' column")
  }
  
  # Standardize column names
  if ("chromosome" %in% names(df)) {
    df$chr <- df$chromosome
  }
  
  if ("position" %in% names(df) && !"end" %in% names(df)) {
    df$end <- df$position
  } else if (!"end" %in% names(df)) {
    df$end <- df$start
  }
  
  # Create GRanges
  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end)
  )
  
  return(gr)
}

#' Get TxDb object for genome
#'
#' @param genome Genome version (e.g., "hg38", "hg19", "mm10")
#' @return TxDb object
#' @export
get_txdb <- function(genome = "hg38") {
  
  load_annotation_libraries()
  
  # Map genome to TxDb package
  txdb_map <- list(
    "hg38" = "TxDb.Hsapiens.UCSC.hg38.refGene",
    "hg19" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "mm10" = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "mm9" = "TxDb.Mmusculus.UCSC.mm9.knownGene"
  )
  
  txdb_name <- txdb_map[[genome]]
  
  if (is.null(txdb_name)) {
    stop("Unsupported genome: ", genome, 
         "\nSupported genomes: ", paste(names(txdb_map), collapse = ", "))
  }
  
  # Load TxDb
  if (!requireNamespace(txdb_name, quietly = TRUE)) {
    message("Installing TxDb package: ", txdb_name)
    BiocManager::install(txdb_name, ask = FALSE)
  }
  
  txdb <- get(txdb_name)
  
  message("Loaded TxDb: ", txdb_name)
  
  return(txdb)
}

#' Run complete annotation and enrichment
#'
#' @param dmr_data Data frame with DMR/DMS results
#' @param config Configuration list
#' @return List with annotation and enrichment results
#' @export
run_annotation <- function(dmr_data, config) {
  
  message("\n")
  message("=" <- rep("=", 60))
  message("Annotation and Enrichment")
  message("=" <- rep("=", 60))
  
  # Get TxDb
  txdb <- get_txdb(config$genome)
  
  # Get organism
  organism <- if (grepl("^hg|^GRCh", config$genome, ignore.case = TRUE)) {
    "hsa"
  } else if (grepl("^mm|^GRCm", config$genome, ignore.case = TRUE)) {
    "mmu"
  } else {
    stop("Cannot determine organism for genome: ", config$genome)
  }
  
  # Get org database name
  org_db <- switch(organism,
                   "hsa" = "org.Hs.eg.db",
                   "mmu" = "org.Mm.eg.db")
  
  # Load org database if available
  if (!requireNamespace(org_db, quietly = TRUE)) {
    message("Installing organism database: ", org_db)
    BiocManager::install(org_db, ask = FALSE)
  }
  
  results_dir <- config$output_dirs$results
  plots_dir <- file.path(results_dir, "plots")
  
  # Annotate DMRs
  message("\nAnnotating DMRs...")
  annotated <- annotate_with_chipseeker(
    dmr_data,
    txdb = txdb,
    annoDb = get(org_db)
  )
  
  # Save annotated data
  annotated_file <- file.path(results_dir, "DMR_annotated.csv")
  write.csv(annotated, annotated_file, row.names = FALSE)
  message("Saved annotated DMRs:", annotated_file)
  
  # Extract genes for enrichment
  # Filter by q-value and methylation difference
  significant_dmrs <- annotated[annotated$pvalue < 0.05 & 
                                 annotated$meth.diff < -10, ]
  genes <- unique(significant_dmrs$geneId)
  
  message("\nGenes for enrichment:", length(genes))
  
  # Run GO enrichment
  go_result <- run_go_enrichment(genes, organism = organism, ontology = "ALL")
  
  if (!is.null(go_result)) {
    # Save GO results
    go_file <- file.path(results_dir, "GO_enrichment.csv")
    write.csv(go_result@result, go_file, row.names = FALSE)
    message("Saved GO enrichment results:", go_file)
    
    # Plot GO
    go_plot <- file.path(plots_dir, "GO_enrichment.png")
    plot_enrichment(go_result, plot_type = "dot", output_file = go_plot)
  }
  
  # Run KEGG enrichment
  kegg_result <- run_kegg_enrichment(genes, organism = organism)
  
  if (!is.null(kegg_result)) {
    # Save KEGG results
    kegg_file <- file.path(results_dir, "KEGG_enrichment.csv")
    write.csv(kegg_result@result, kegg_file, row.names = FALSE)
    message("Saved KEGG enrichment results:", kegg_file)
    
    # Plot KEGG
    kegg_plot <- file.path(plots_dir, "KEGG_enrichment.png")
    plot_enrichment(kegg_result, plot_type = "dot", output_file = kegg_plot)
  }
  
  # Return results
  results <- list(
    annotated = annotated,
    go_result = go_result,
    kegg_result = kegg_result
  )
  
  return(results)
}
