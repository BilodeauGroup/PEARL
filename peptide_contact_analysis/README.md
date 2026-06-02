<contents of README.md>
# 🧬 Peptide Contact Matrix Analysis

This repository computes and visualizes **residue-residue contact matrices** for pentapeptides from molecular dynamics trajectories.

## 🧠 Overview

- **compute_contact_matrix.py** → Uses MDAnalysis to compute contact frequency matrices between peptide residues.
- **plot_contact_heatmaps.py** → Reads CSV outputs and generates publication-quality heatmaps.

## 🗂️ Directory Layout
peptide_contact_analysis/
├── data/contact_matrices/ # CSV outputs
├── results/heatmaps/ # PNG heatmaps
├── scripts/
│ ├── compute_contact_matrix.py
│ └── plot_contact_heatmaps.py
├── requirements.txt
└── README.md
```bash
git clone https://github.com/YOUR_USERNAME/peptide_contact_analysis.git
cd peptide_contact_analysis
pip install -r requirements.txt
🚀 Usage
1️⃣ Compute Contact Matrices
python scripts/compute_contact_matrix.py

2️⃣ Plot Heatmaps
python scripts/plot_contact_heatmaps.py

📊 Output Example

Each sequence–replicate combination generates:

data/contact_matrices/{seq}/{rep}_contact_matrix.csv

results/heatmaps/{seq}_{rep}_contact_heatmap.png
