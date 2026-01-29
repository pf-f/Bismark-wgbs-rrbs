# Installation Guide

This guide will help you set up the Methylation Analysis Toolkit for WGBS and RRBS analysis.

## Prerequisites

### System Requirements
- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+) or macOS 10.15+
- **RAM**: Minimum 16GB, recommended 32GB+ for WGBS
- **CPU**: Multi-core processor (4+ cores recommended)
- **Storage**: 100GB+ free space (depending on data size)

### Software Requirements
- **Conda**: Miniconda or Anaconda (Python 3.8+)
- **R**: R 4.0+ (with Bioconductor)
- **Java**: OpenJDK 8+ (for Bismark)

## Step-by-Step Installation

### 1. Install Conda (if not already installed)

**Linux:**
```bash
# Download Miniconda
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py310_24.1.2-0-Linux-x86_64.sh

# Install Miniconda
bash Miniconda3-py310_24.1.2-0-Linux-x86_64.sh

# Follow the prompts (accept license, choose install location)

# Activate conda
source ~/.bashrc
```

**macOS:**
```bash
# Download Miniconda
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py310_24.1.2-0-MacOSX-x86_64.sh

# Install Miniconda
bash Miniconda3-py310_24.1.2-0-MacOSX-x86_64.sh
```

### 2. Clone the Repository

```bash
git clone https://github.com/yourusername/methylation_analysis.git
cd methylation_analysis
```

### 3. Create Conda Environment

**For WGBS:**
```bash
conda env create -f environments/wgbs_env.yml
conda activate methylation_wgbs
```

**For RRBS:**
```bash
conda env create -f environments/rrbs_env.yml
conda activate methylation_rrbs
```

This will install:
- Bismark
- Bowtie2
- Samtools
- FastQC
- MultiQC
- Trimmomatic (WGBS) or Trim Galore (RRBS)
- Python 3.10
- R 4.3

### 4. Install R Dependencies

```bash
Rscript install.R
```

This will install all required R packages including:
- **Bioconductor**: methylKit, genomation, ChIPseeker, clusterProfiler
- **CRAN**: ggplot2, dplyr, tidyr, optparse, yaml
- **Databases**: org.Hs.eg.db, org.Mm.eg.db, TxDb packages
- **Optional**: RIdeogram, gggenes, ggrepel

**Note**: Installing TxDb packages may take 10-20 minutes.

### 5. Install Python Dependencies

```bash
pip install -r requirements.txt
```

### 6. Verify Installation

**Check conda environment:**
```bash
conda list  # Verify installed packages
```

**Check R packages:**
```bash
R -e "library(methylKit); library(genomation); library(clusterProfiler)"
```

**Check command-line tools:**
```bash
which bismark
which samtools
which fastqc
```

**Run Python environment check:**
```bash
python src/methylation_analysis/__init__.py
```

## Download Reference Genome

### Human Genome (hg38)

```bash
# Create genome directory
mkdir -p data/ref/GRCh38

# Download genome FASTA
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz \
  -O data/ref/GRCh38/GRCh38.fa.gz

# Unzip
gunzip data/ref/GRCh38/GRCh38.fa.gz

# Rename
mv data/ref/GRCh38/GRCh38.fa data/ref/GRCh38/GRCh38.fa
```

### Mouse Genome (mm10)

```bash
# Create genome directory
mkdir -p data/ref/mm10

# Download genome FASTA
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz \
  -O data/ref/mm10/mm10.fa.gz

# Unzip
gunzip data/ref/mm10/mm10.fa.gz
```

### Build Bismark Genome Indices

```bash
# Activate environment
conda activate methylation_wgbs  # or methylation_rrbs

# Build genome indices (this may take several hours)
bismark_genome_preparation --verbose data/ref/GRCh38/

# Or for mouse
bismark_genome_preparation --verbose data/ref/mm10/
```

## Download Annotation Files

### GTF/Annotation Files

```bash
# Create bed directory
mkdir -p data/ref/bed

# Download human annotation
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Hg38_RefSeq_BED.tar.gz \
  -O data/ref/bed/hg38_annotation.tar.gz

# Extract
tar -xzf data/ref/bed/hg38_annotation.tar.gz -C data/ref/bed/
```

## Running Tests

**Test WGBS setup:**
```bash
Rscript tests/test_wgbs.R
```

**Test RRBS setup:**
```bash
Rscript tests/test_rrbs.R
```

## Troubleshooting

### Issue: Conda installation fails

**Solution:**
```bash
# Update conda
conda update conda

# Try using different mirror
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
```

### Issue: R package installation fails

**Solution:**
```bash
# Update Bioconductor
R
> BiocManager::install(update = TRUE)

# Install specific package
BiocManager::install("methylKit", ask = FALSE, force = TRUE)
```

### Issue: Out of memory during alignment

**Solution:**
- Reduce number of threads in config file
- Increase system swap space
- Use smaller genome (e.g., hg19 instead of hg38)
- Reduce buffer_size in config

### Issue: Bismark not found after conda activate

**Solution:**
```bash
# Make sure to activate the environment
conda activate methylation_wgbs

# Check PATH
which bismark

# If not found, add to PATH manually
export PATH=$CONDA_PREFIX/bin:$PATH
```

### Issue: Genome preparation takes too long

**Solution:**
```bash
# Use more threads (if available)
bismark_genome_preparation --parallel 8 --verbose data/ref/GRCh38/

# Or use pre-built indices if available
```

## Environment Variables (Optional)

You can set these in your `~/.bashrc` or `~/.zshrc`:

```bash
# Default threads
export METHYLATION_THREADS=8

# Default memory (GB)
export METHYLATION_MEMORY=20

# Default output directory
export METHYLATION_OUTPUT=outputs
```

Then use in scripts:
```bash
Rscript scripts/run_analysis.R --threads $METHYLATION_THREADS
```

## Next Steps

Once installation is complete:

1. **Prepare your data**: Place FASTQ files in `data/raw_fastq/`
2. **Create sample info**: Create `sample_info.csv` with sample information
3. **Run analysis**: Execute the pipeline using `run_analysis.R`
4. **Review results**: Check `outputs/05_results/` for final results

For detailed usage instructions, see [README.md](README.md).

## Support

If you encounter issues:
1. Check the log files in `outputs/log/`
2. Review error messages carefully
3. Consult the main [README.md](README.md)
4. Open an issue on GitHub with:
   - Your operating system
   - Conda environment info (`conda list`)
   - R version (`R --version`)
   - Error message
   - Steps to reproduce
