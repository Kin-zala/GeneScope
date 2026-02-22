# 🧬 GeneScope

**GeneScope** is a Streamlit app for interactive gene expression visualization and exploratory analysis.

It supports:  
- Uploading gene expression CSVs  
- Standard visualization (boxplot & heatmap)  
- PCA for large datasets with KMeans clustering  
- Log2 transformation option  
- Downloadable plots and clustering results  

---

## 🚀 Quick Start

1. Clone the repo:

```bash
git clone https://github.com/Kin-zala/GeneScope.git
cd GeneScope 
```
2. Install dependencies:
```bash
pip install -r requirements.txt
```
3. Run the app
```bash
streamlit run genescope_app.py
```

## 📂 Input Format

First column: Gene names

Remaining columns: Numeric expression values

Rows: Genes, Columns: Samples

#### Example CSV:

| Gene | Sample_1 | Sample_2 | Sample_3 |
|:-----|:--------:|:--------:|:--------:|
| TP53 |  120 |  98  |  101 |
| EGFR |  300 |  287 |  310 |

## 🧪 Demo Dataset

Click "Load Demo Dataset" in the app to try a random dataset with 2000 genes × 50 samples.

## ⚡ Features

Top variable genes selection

Boxplot & heatmap visualization

PCA + KMeans clustering

Interactive Plotly charts

Downloadable plots (PNG) and tables (CSV)