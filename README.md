## PCA Analysis on Breast Cancer Gene Expression Data

A clean and beginner-friendly Python project for visualizing gene expression patterns and exploring dimensionality reduction using **Principal Component Analysis (PCA)** on a breast cancer dataset.
This project recreates important biological insights using simple data science techniques and clear visualizations. 

---

## Project Overview

This project works on a breast cancer gene expression dataset and focuses on:

* Comparing the expression of **XBP1** and **GATA3** genes
* Visualizing the separation between **ER+** and **ER−** cancer samples
* Applying **PCA (Principal Component Analysis)** to reduce high-dimensional gene data into simpler visual representations
* Generating publication-style plots for analysis

The implementation is written completely in Python using common data science libraries like NumPy, Pandas, Matplotlib, and Scikit-learn.

---

## Features

✅ Loads compressed gene expression datasets
✅ Cleans and preprocesses data
✅ Extracts important gene markers (XBP1 & GATA3)
✅ Creates scatter plots for biological comparison
✅ Performs PCA on the full dataset
✅ Visualizes PCA projections with class separation
✅ Saves high-quality figures automatically

---

## Technologies Used

* Python 3
* NumPy
* Pandas
* Matplotlib
* Scikit-learn

---

## Dataset Information

The project expects two input files:

| File              | Description                         |
| ----------------- | ----------------------------------- |
| `filtered_tsv.gz` | Gene expression matrix              |
| `class.tsv`       | ER status labels (1 = ER+, 0 = ER−) |

---

##  How to Run

### 1. Clone the repository

```bash
git clone <your-repository-link>
cd <repository-folder>
```

### 2. Install dependencies

```bash
pip install numpy pandas matplotlib scikit-learn
```

### 3. Place dataset files

Put the following files inside the required input directory:

* `filtered_tsv.gz`
* `class.tsv`

### 4. Run the script

```bash
python pca_analysis.py
```

---

## Output

The script generates two figures:

### Figure 1a — XBP1 vs GATA3 Scatter Plot

Shows how ER+ and ER− samples differ based on expression levels of two important genes.

### Figure 1c — PCA Projection Plot

Displays samples projected onto the first principal component (PC1), helping visualize clustering and variance separation.

Generated images:

```bash
figure_1a_scatter.png
figure_1c_pca.png
```

---

## Understanding the Workflow

### Step 1 — Load Data

The gene expression matrix and class labels are loaded using Pandas.

### Step 2 — Extract Important Genes

The script selects:

* **XBP1** → Gene ID `4404`
* **GATA3** → Gene ID `4359`

### Step 3 — Visualization

Scatter plots are created to compare gene expression patterns between cancer subtypes.

### Step 4 — Standardization

Data is standardized before PCA because genes may have different expression scales.

### Step 5 — PCA

PCA reduces thousands of gene dimensions into principal components while preserving maximum variance.

### Step 6 — Plot PCA Results

Samples are projected onto PC1 to observe clustering between ER+ and ER− groups.

---

## Learning Outcomes

This project is useful for understanding:

* Gene expression analysis
* Data preprocessing
* PCA and dimensionality reduction
* Biological data visualization
* Classification pattern exploration

---


