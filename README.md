# DNA Methylation Analysis Toolkit

A comprehensive toolkit for DNA methylation analysis supporting both **WGBS** (Whole Genome Bisulfite Sequencing) and **RRBS** (Reduced Representation Bisulfite Sequencing) workflows.


### WGBS Workflow
```
FastQ → QC → Trimming → Bismark Alignment → Deduplication → 
Methylation Extraction → Differential Analysis → Annotation → Visualization
```

### RRBS Workflow
```
FastQ → QC → Trim Galore → Bismark Alignment → 
Methylation Extraction → Differential Analysis → Annotation → Visualization
```

## Command Line Options

```
Usage: run_analysis.R [options]

Options:
  -t, --type             Analysis type: WGBS or RRBS [default: WGBS]
  -g, --genome           Genome version (e.g., hg38, mm10) [default: hg38]
  -i, --input           Input directory with FASTQ files [default: data/raw_fastq]
  -c, --config          Configuration file (YAML)
  -o, --output          Output directory [default: outputs]
  -s, --steps           Comma-separated steps: qc,align,dedup,extract,diff,annotate,viz
  --threads              Number of threads [default: 8]
  --sample-info          Sample information file (CSV)
  --min-cov             Minimum coverage threshold [default: 10]
  --diff-thresh         Methylation difference threshold (%) [default: 25]
  --qvalue              Q-value threshold [default: 0.01]
  -v, --verbose         Print verbose output [default: TRUE]
  --version             Print version and exit
```

## Project Structure

```
methylation_analysis/
├── src/
│   └── methylation_analysis/
│       ├── __init__.py              # Python package initialization
│       └── R/
│           ├── config.R             # Configuration management
│           ├── qc_processor.R       # QC and preprocessing
│           ├── aligner.R            # Bismark alignment
│           ├── deduplication.R      # Deduplication (WGBS only)
│           ├── methylation_extractor.R # Methylation extraction
│           ├── analyzer.R           # Differential analysis
│           ├── annotator.R          # Annotation and enrichment
│           ├── visualizer.R         # Visualization
│           └── pipeline.R           # Main pipeline controller
├── scripts/
│   └── run_analysis.R               # Command-line entry point
├── examples/
│   ├── WGBS_tutorial.Rmd            # Complete WGBS tutorial
│   ├── RRBS_tutorial.Rmd            # Complete RRBS tutorial
├── tests/
│   ├── test_wgbs.R                 # WGBS unit tests
│   └── test_rrbs.R                 # RRBS unit tests
├── data/
│   ├── raw_fastq/                   # Input FASTQ files
│   ├── ref/                        # Reference genomes
│   │   ├── GRCh38/                 # Human genome
│   │   └── mm10/                   # Mouse genome
│   └── bed/                        # Annotation BED files
├── outputs/                        # Output directory
│   ├── 01_clip_reads/             # Trimmed FASTQ files
│   ├── 02_meth_bam/              # Bismark BAM files
│   ├── 03_dedu_meth/             # Deduplicated BAMs (WGBS)
│   ├── 04_methinfo/              # Methylation information
│   └── 05_results/               # Final results
│       ├── methylkit/              # methylKit analysis
│       ├── plots/                  # Visualization plots
│       └── *.csv                  # Annotated results
├── environments/
│   ├── wgbs_env.yml              # WGBS conda environment
│   └── rrbs_env.yml              # RRBS conda environment
├── requirements.txt                 # Python dependencies
├── install.R                       # R dependencies installation
├── pyproject.toml                  # Python package config
└── README.md
```

## Analysis Steps

### 1. Quality Control
- **FastQC**: Raw read quality assessment
- **Trimmomatic/Trim Galore**: Adapter trimming and quality filtering
- **MultiQC**: Aggregate QC report

### 2. Alignment
- **Bismark**: Bisulfite-seq aligner using Bowtie2
- **Genome preparation**: Build Bismark genome indices
- **Parallel processing**: Multi-threaded alignment

### 3. Deduplication (WGBS only)
- **deduplicate_bismark**: Remove PCR duplicates
- **Note**: Skipped for RRBS data

### 4. Methylation Extraction
- **bismark_methylation_extractor**: Extract CpG, CHG, CHH methylation
- **Coverage files**: Generate per-context coverage files
- **Cytosine report**: Comprehensive methylation statistics

### 5. Differential Analysis
- **methylKit**: DMR/DMS calculation
- **Coverage filtering**: Remove low-coverage sites
- **Statistical tests**: Fisher's exact test / logistic regression
- **Clustering**: Sample correlation and PCA

### 6. Annotation
- **genomation**: Gene region annotation
- **ChIPseeker**: Peak annotation with gene information
- **TxDb**: Gene model annotations

### 7. Enrichment Analysis
- **clusterProfiler**: GO and KEGG enrichment
- **Biological interpretation**: Functional analysis of DMR-associated genes

### 8. Visualization
- **Genome-wide plots**: Chromosome ideograms
- **Volcano plots**: DMR significance
- **Heatmaps**: Sample methylation patterns
- **Enrichment plots**: GO/KEGG term visualization

## Version History

