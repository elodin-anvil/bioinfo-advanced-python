# Bioinformatics Advanced Python Project - Instructions

## Introduction
Welcome to your guide for the Bioinformatics Advanced Python project! This document will walk you through the theory behind single-cell RNA sequencing (scRNA-seq) analysis, your specific role in the project, and how to set up and run the code provided.

We are building a collaborative pipeline to process raw gene expression data into meaningful biological insights (clusters of cell types).

---

## 1. The Big Picture: scRNA-seq Pipeline

Imagine you have a spreadsheet where:
*   **Rows** are different genes (thousands of them).
*   **Columns** are individual cells (thousands of them).
*   The **values** represent how active each gene is in each cell (expression count).

Our goal is to group these cells based on their gene activity patterns. Cells with similar patterns are likely the same cell type (e.g., T-cells, B-cells).

### The Pipeline Steps:
1.  **Data Loading & Preprocessing (You & Simon):** Import the messy raw data, clean it up, and get it ready for analysis.
2.  **Dimensionality Reduction (Stephan):** Simplify the data (PCA) so computers can handle it.
3.  **Clustering (Raeez, Stephan, Seb):** Group similar cells together.
4.  **Visualization (Samishka):** Create plots (UMAP) to see these groups.

---

## 2. Your Role: Data Ingestion & Preprocessing

Your job is the foundation of the entire project. If the data isn't loaded and cleaned correctly, everything downstream will fail or give wrong results.

### Key Tasks:
1.  **Load the Data:** Read the CSV files.
2.  **Transpose:** Ensure the data is in the format `(Cells x Genes)`.
    *   *Why?* Most raw data comes as `(Genes x Cells)`, but analysis tools like `Scanpy` expect `(Cells x Genes)`.
3.  **Quality Control (QC):**
    *   Remove bad cells (e.g., empty or dead cells).
    *   Remove uninformative genes (genes that aren't expressed anywhere).
4.  **Normalization:**
    *   Scale the data so that cells with more total reads don't overshadow cells with fewer reads (just because of sequencing depth).
    *   We use `log1p` (logarithm + 1) to handle the wide range of expression values.

---

## 3. Getting Started

### Prerequisites
You need Python installed. We will use a virtual environment to manage dependencies.

### Step 1: Clone the Repository
Open your terminal (Command Prompt or PowerShell on Windows, Terminal on Mac/Linux) and run:

```bash
git clone https://github.com/elodin-anvil/bioinfo-advanced-python.git
cd bioinfo-advanced-python
```

### Step 2: Set Up Environment
Create a virtual environment (keeps your project libraries separate from your system):

**Windows:**
```bash
python -m venv .venv
.venv\Scripts\activate
```

**Mac/Linux:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### Step 3: Install Dependencies
We use `scanpy`, `pandas`, and other libraries. Install them with:

```bash
pip install -r requirements.txt
```

---

## 4. The Code: `src/preprocessing/data_pipeline.py`

I have created a script for you in `src/preprocessing/data_pipeline.py`. Let's look at what it does.

### The `load_and_preprocess_data` function:
This is the core function you will use.

```python
def load_and_preprocess_data(file_path, min_genes=200, min_cells=3):
    # ... code ...
```

*   **`file_path`**: The location of your raw CSV file.
*   **`sc.read_csv(file_path).T`**: Reads the CSV and Transposes it (.T).
*   **`sc.pp.normalize_total`**: Normalizes the data.
*   **`sc.pp.log1p`**: Log-transforms the data.
*   **`sc.pp.highly_variable_genes`**: Identifies genes that vary the most (these are the most useful for clustering).

### Running the Script
You can test the script right now! It includes a small test block at the bottom that creates dummy data if no real data is found.

Run it from the root directory:
```bash
python src/preprocessing/data_pipeline.py
```

You should see output like:
```
Created dummy data at data/raw/test_data.csv
File validated successfully: data/raw/test_data.csv
Initial data shape (Cells x Genes): (5, 10)
...
Data preprocessing complete.
```

---

## 5. Next Steps for You

1.  **Get the Real Data:** Ask your teammates (or check the shared drive) for the actual `.csv` files.
2.  **Place Data:** Put them in the `data/raw/` folder (create it if it's missing).
3.  **Update the Script:** Modify the `if __name__ == "__main__":` block in `src/preprocessing/data_pipeline.py` to point to the real file name instead of `test_data.csv`.
4.  **Run & Verify:** Run the script again. If it prints "Data preprocessing complete" and the shape looks correct (e.g., thousands of cells and genes), you are done!
5.  **Commit & Push:**
    ```bash
    git add src/preprocessing/data_pipeline.py
    git commit -m "Added data preprocessing pipeline"
    git push origin main
    ```

Good luck! Let me know if you have questions about the "why" or "how" of any step.
