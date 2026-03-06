# Changelog

All notable changes to microbiome-disease-pipeline will be documented in this file.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Kraken2 classification and Bracken abundance estimation modules
- Alpha diversity metrics (Shannon, Simpson, Chao1, observed richness)
- Beta diversity (Bray-Curtis, Jaccard) with PCoA ordination
- Unit tests for diversity calculations
- CI/CD pipeline (GitHub Actions)
- Project documentation and README with architecture diagram

### Planned
- MetaPhlAn4 taxonomic profiling integration
- HUMAnN3 functional annotation pipeline
- KneadData host DNA removal
- Differential abundance testing
- ML disease prediction models
- Visualization module
- Snakemake workflow for reproducibility
