"""
Methylation Analysis Toolkit

A comprehensive toolkit for DNA methylation analysis supporting both 
WGBS and RRBS workflows.

This package provides:
- WGBS (Whole Genome Bisulfite Sequencing) analysis
- RRBS (Reduced Representation Bisulfite Sequencing) analysis
- Quality control and preprocessing
- Alignment with Bismark
- Differential methylation analysis
- Annotation and functional enrichment
- Visualization

Installation:
    # Create conda environment
    conda env create -f environments/wgbs_env.yml  # For WGBS
    conda env create -f environments/rrbs_env.yml  # For RRBS
    
    # Install R dependencies
    Rscript install.R
    
    # Install Python dependencies
    pip install -r requirements.txt

Usage:
    # Run WGBS analysis
    Rscript scripts/run_analysis.R -t WGBS -i data/raw_fastq -o outputs/
    
    # Run RRBS analysis
    Rscript scripts/run_analysis.R -t RRBS -i data/raw_fastq -o outputs/

For detailed tutorials, see examples/WGBS_tutorial.Rmd and examples/RRBS_tutorial.Rmd
"""

__version__ = "1.0.0"
__author__ = "Methylation Analysis Team"

# Python utility functions (for future extensions)

def get_config_path():
    """
    Get the default configuration file path.
    
    Returns:
        str: Path to default configuration file
    """
    import os
    return os.path.join(os.path.dirname(__file__), "config.yaml")


def validate_environment():
    """
    Validate that required tools are installed and accessible.
    
    Returns:
        dict: Dictionary with tool names and installation status
    """
    import subprocess
    
    required_tools = {
        'bismark': 'bismark --version',
        'bismark_genome_preparation': 'bismark_genome_preparation --version',
        'deduplicate_bismark': 'deduplicate_bismark --version',
        'bismark_methylation_extractor': 'bismark_methylation_extractor --version',
        'samtools': 'samtools --version',
        'fastqc': 'fastqc --version',
    }
    
    results = {}
    
    for tool, cmd in required_tools.items():
        try:
            result = subprocess.run(cmd, shell=True, 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=5)
            results[tool] = {
                'installed': result.returncode == 0,
                'version': result.stdout.strip() if result.returncode == 0 else None
            }
        except subprocess.TimeoutExpired:
            results[tool] = {'installed': False, 'version': None}
        except Exception as e:
            results[tool] = {'installed': False, 'version': None, 'error': str(e)}
    
    return results


def print_environment_status():
    """Print the status of required tools."""
    status = validate_environment()
    
    print("=" * 60)
    print("Environment Status")
    print("=" * 60)
    
    for tool, info in status.items():
        if info['installed']:
            print(f"  {tool:30} ✓   {info['version']}")
        else:
            print(f"  {tool:30} ✗   Not installed")
    
    print("=" * 60)


if __name__ == "__main__":
    print_environment_status()
