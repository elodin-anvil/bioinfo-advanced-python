# Bioinformatics Advanced Python Project
## scRNA-seq Pipeline Walkthrough

Your Guide to the Galaxy of Cells 🌌

---

# The Goal 🎯

We want to group cells based on their gene activity.

- **Input:** Raw gene expression data (count matrix).
- **Process:** Clean -> Normalize -> Reduce Dimensions -> Cluster.
- **Output:** Clusters of cell types (e.g., T-cells vs B-cells).

---

# Your Role: The Foundation 🏗️

**Data Ingestion & Preprocessing**

If you get this wrong, the rest of the team fails.
Your job is to feed clean, normalized data to the clustering algorithms.

---

# Step 1: Loading & Transposing 🔄

**Problem:**
Raw data often comes as `Genes x Cells` (genes are rows).
Python tools (`Scanpy`) expect `Cells x Genes` (cells are rows).

**Solution:**
We transpose the matrix immediately after loading.

```python
adata = sc.read_csv(file_path).T
```

---

# Step 2: Quality Control (QC) 🧹

**Why?**
- **Dead cells:** Have very few expressed genes.
- **Empty droplets:** The sequencer read a droplet with no cell inside.
- **Doublets:** Two cells stuck together (looks like one super-cell).

**Action:**
Filter out cells with too few genes (`min_genes=200`).
Filter out genes seen in too few cells (`min_cells=3`).

---

# Step 3: Normalization ⚖️

**Why?**
Sequencing depth varies.
- Cell A: 10,000 reads
- Cell B: 5,000 reads

If Cell A has 2x more reads for Gene X, is it biologically more active?
Maybe not. It might just be sequenced deeper.

**Action:**
Normalize every cell to have the same total count (e.g., 10,000).

```python
sc.pp.normalize_total(adata, target_sum=1e4)
```

---

# Step 4: Log Transformation 📉

**Why?**
Gene expression follows a power law.
Some genes are expressed 1000x more than others.
This huge range dominates statistical analysis.

**Action:**
Take the log (`log1p` = log(x+1)) to bring values into a comparable range.

```python
sc.pp.log1p(adata)
```

---

# Step 5: Feature Selection 🌟

**Why?**
We have 20,000+ genes.
Most are "housekeeping" genes (always on) or noise (rarely on).
They don't help distinguish cell types.

**Action:**
Select "Highly Variable Genes" (HVGs) — genes that differ strongly between cells.
These hold the biological signal.

```python
sc.pp.highly_variable_genes(adata)
```

---

# Your Pipeline Code 💻

I've written `src/preprocessing/data_pipeline.py` for you.

**How to run it:**
1.  Place raw `.csv` files in `data/raw/`.
2.  Update the filename in the script.
3.  Run: `python src/preprocessing/data_pipeline.py`

---

# Next Steps 🚀

1.  **Clone the Repo:** `git clone ...`
2.  **Install:** `pip install -r requirements.txt`
3.  **Run:** Execute the pipeline script.
4.  **Verify:** Check the output shape (should be `Cells x Genes`).

Good luck! You're building the bedrock of this analysis.
