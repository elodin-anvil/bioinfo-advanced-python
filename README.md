# Bioinformatics Advanced Python Project

## Overview
This repository contains the codebase for our group project on single-cell RNA sequencing (scRNA-seq) analysis. The pipeline is designed to ingest raw gene expression data, preprocess it, perform dimensionality reduction and clustering, and visualize the results.

## Project Structure
```
bioinfo-advanced-python/
├── data/
│   ├── raw/            # Place raw .csv files here
│   └── processed/      # Cleaned and normalized data
├── src/
│   ├── preprocessing/  # Data loading and QC scripts
│   ├── clustering/     # Clustering algorithms (Leiden, DBSCAN, K-Means)
│   └── visualization/  # Plotting functions (UMAP)
├── tests/              # Unit tests
├── INSTRUCTIONS.md     # Detailed guide for team members
└── requirements.txt    # Python dependencies
```

## Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/elodin-anvil/bioinfo-advanced-python.git
    cd bioinfo-advanced-python
    ```

2.  **Create a virtual environment:**
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Data Preprocessing
To load and preprocess the raw data:
```bash
python src/preprocessing/data_pipeline.py
```
This script handles:
*   File validation
*   Data transposition (Cells x Genes)
*   QC (filtering low-quality cells/genes)
*   Normalization (log1p)

## Team
*   **Data Prep:** Jonny & Simon
*   **Clustering:** Raeez, Stephan, Seb
*   **Visualization:** Samishka
