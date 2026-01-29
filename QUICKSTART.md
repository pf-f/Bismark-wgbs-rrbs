# Quick Start Guide

This guide will help you run your first methylation analysis with minimal setup.

## 1. Minimal Setup (5 minutes)

### Install environment (choose one)

**WGBS:**
```bash
conda env create -f environments/wgbs_env.yml
conda activate methylation_wgbs
```

**RRBS:**
```bash
conda env create -f environments/rrbs_env.yml
conda activate methylation_rrbs
```

### Install R packages (2 minutes)
```bash
Rscript install.R
```

## 2. Prepare Data (1 minute)

### Create sample info file

Create `sample_info.csv` with your sample information:

```csv
run_accession,sample_alias,treatment,treatment_code
SRR001,sample1,control,0
SRR002,sample2,control,0
SRR003,sample3,treated,1
SRR004,sample4,treated,1
```

### Place FASTQ files

```bash
# Create data directory
mkdir -p data/raw_fastq

# Copy your FASTQ files
cp /path/to/your/SRR001_1.fastq.gz data/raw_fastq/
cp /path/to/your/SRR001_2.fastq.gz data/raw_fastq/
# ... repeat for all samples
```

## 3. Run Analysis (2 minutes)

### WGBS analysis

```bash
Rscript scripts/run_analysis.R \
  -t WGBS \
  -i data/raw_fastq \
  -o outputs \
  --sample-info sample_info.csv
```

### RRBS analysis

```bash
Rscript scripts/run_analysis.R \
  -t RRBS \
  -i data/raw_fastq \
  -o outputs \
  --sample-info sample_info.csv
```

### Run specific steps only

If you want to run only certain steps:

```bash
# Only run QC and alignment
Rscript scripts/run_analysis.R \
  -t WGBS \
  -i data/raw_fastq \
  -o outputs \
  -s qc,align
```

Available steps: `qc,align,dedup,extract,diff,annotate,viz`

## 4. Check Results

### Quality control
```bash
# View QC reports
ls outputs/fastqc/
ls outputs/multiqc_report.html  # Open in browser
```

### Alignment results
```bash
# Check alignment reports
ls outputs/02_meth_bam/*_report.txt
cat outputs/02_meth_bam/*_report.txt
```

### Differential methylation
```bash
# View DMR results
ls outputs/05_results/methylkit/
head outputs/05_results/methylkit/hyper_dms.csv
head outputs/05_results/methylkit/hypo_dms.csv
```

### Annotation results
```bash
# View annotated DMRs
head outputs/05_results/DMR_annotated.csv

# View enrichment results
head outputs/05_results/GO_enrichment.csv
head outputs/05_results/KEGG_enrichment.csv
```

### Plots
```bash
# View all plots
ls outputs/05_results/plots/
```

## 5. Common Use Cases

### Analyze only a subset of samples

Edit `sample_info.csv` to include only the samples you want to analyze.

### Compare specific groups

Edit `sample_info.csv` to set appropriate treatment codes:
- `0` = Control group
- `1` = Treatment group 1
- `2` = Treatment group 2

### Adjust sensitivity

Change thresholds in config file or command line:

```bash
# More sensitive analysis (more DMRs, higher false positive rate)
Rscript scripts/run_analysis.R \
  -t WGBS \
  --diff-thresh 15 \
  --qvalue 0.05

# More conservative analysis (fewer DMRs, lower false positive rate)
Rscript scripts/run_analysis.R \
  -t WGBS \
  --diff-thresh 30 \
  --qvalue 0.001
```

### Use custom configuration

1. Copy config file:
```bash
cp config_wgbs.yaml my_config.yaml
```

2. Edit parameters in `my_config.yaml`

3. Run with custom config:
```bash
Rscript scripts/run_analysis.R -c my_config.yaml
```

## Troubleshooting

### "No FASTQ files found"

Make sure:
- FASTQ files are in `data/raw_fastq/`
- Files have `.fastq.gz` extension
- Path is correct

### "sample_info is required for differential analysis"

Create `sample_info.csv` with required columns.

### "Bismark not found"

Activate conda environment:
```bash
conda activate methylation_wgbs  # or methylation_rrbs
```

### "Out of memory"

Reduce number of threads:
```bash
Rscript scripts/run_analysis.R --threads 4
```

### "R package not installed"

Run the installation script again:
```bash
Rscript install.R
```

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Check [INSTALL.md](INSTALL.md) for complete installation guide
- Explore [examples/WGBS_tutorial.Rmd](examples/WGBS_tutorial.Rmd) for detailed tutorial
- Run tests to verify your setup:
  ```bash
  Rscript tests/test_wgbs.R
  Rscript tests/test_rrbs.R
  ```

## Getting Help

If you encounter issues:

1. Check the log file: `outputs/log/pipeline_*.log`
2. Review error messages carefully
3. Consult the main [README.md](README.md)
4. Open an issue on GitHub with your error message

## Example: Complete Workflow

Here's a complete example from start to finish:

```bash
# 1. Setup
git clone https://github.com/yourusername/methylation_analysis.git
cd methylation_analysis
conda env create -f environments/wgbs_env.yml
conda activate methylation_wgbs
Rscript install.R

# 2. Download test data (example)
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/001/SRR12334229/SRR12334229_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/001/SRR12334229/SRR12334229_2.fastq.gz
mv SRR12334229_*.fastq.gz data/raw_fastq/

# 3. Create sample info
cat > sample_info.csv << EOF
run_accession,sample_alias,treatment,treatment_code
SRR12334229,sample1,control,0
EOF

# 4. Run analysis
Rscript scripts/run_analysis.R \
  -t WGBS \
  -i data/raw_fastq \
  -o outputs \
  --sample-info sample_info.csv

# 5. View results
ls outputs/05_results/plots/
```

This will give you:
- Quality control reports
- Aligned BAM files
- Methylation information
- Differential methylation regions
- Annotation with gene information
- Enrichment analysis results
- Visualization plots
