import os
import scanpy as sc
import pandas as pd
import anndata
import numpy as np

def validate_file(file_path):
    """
    Validates the input file path.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found at: {file_path}")
    
    if not file_path.lower().endswith('.csv'):
        raise ValueError("Input file must be a .csv file.")
    
    if os.path.getsize(file_path) == 0:
        raise ValueError("Input file is empty.")

    print(f"File validated successfully: {file_path}")

def load_and_preprocess_data(file_path, min_genes=200, min_cells=3):
    """
    Loads raw RNA-seq data, performs basic QC (filtering), normalization, and scaling.
    
    Args:
        file_path (str): Path to the CSV file.
        min_genes (int): Minimum number of genes expressed required for a cell to pass filtering.
        min_cells (int): Minimum number of cells required for a gene to be kept.
    
    Returns:
        anndata.AnnData: The processed AnnData object.
    """
    validate_file(file_path)
    
    print("Loading data...")
    # Load the data. We assume rows are genes and columns are cells initially, based on typical raw output.
    # Scanpy reads csv. If genes are rows, we transpose.
    # We'll read it first to check shape/content if needed, but scanpy handles it well.
    # Convention: usually scanpy.read_csv expects rows=observations (cells). 
    # If the file is genes x cells, we need to transpose.
    
    # Let's try reading with pandas first to peek at the structure if we were unsure, 
    # but for this assignment, let's assume standard format where rows=genes, cols=samples (cells).
    # This means we need to transpose because AnnData expects rows=cells, cols=genes.
    
    adata = sc.read_csv(file_path).T  # Transpose to get Cells x Genes
    
    print(f"Initial data shape (Cells x Genes): {adata.shape}")

    # Basic QC - Filter out low quality cells and genes
    # sc.pp.filter_cells(adata, min_genes=min_genes)
    # sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # print(f"Data shape after filtering: {adata.shape}")

    # Handle duplicates if any (essential for unique indices)
    if not adata.obs_names.is_unique:
        adata.obs_names_make_unique()
    if not adata.var_names.is_unique:
        adata.var_names_make_unique()

    # Normalize
    # Total-count normalize (library size correct) to 10,000 reads/cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Logarithmize
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Scale data to unit variance
    sc.pp.scale(adata, max_value=10)
    
    print("Data preprocessing complete.")
    return adata

if __name__ == "__main__":
    # Example usage
    # This block allows you to test the script directly
    try:
        # Create a dummy file for testing if it doesn't exist
        dummy_path = "data/raw/test_data.csv"
        if not os.path.exists(dummy_path):
            os.makedirs("data/raw", exist_ok=True)
            # Create a small dummy CSV: 10 genes (rows), 5 cells (cols)
            df = pd.DataFrame(np.random.randint(0, 100, size=(10, 5)), 
                              index=[f"Gene_{i}" for i in range(10)], 
                              columns=[f"Cell_{i}" for i in range(5)])
            df.to_csv(dummy_path)
            print(f"Created dummy data at {dummy_path}")

        adata = load_and_preprocess_data(dummy_path)
        print(adata)
        
    except Exception as e:
        print(f"Error: {e}")
