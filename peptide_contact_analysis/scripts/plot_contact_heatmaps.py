"""
Plot contact matrix heatmaps for all peptide sequences from CSVs.
Outputs figures to ./results/heatmaps/
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm

# --- Configuration ---
sequences = ["KYYYL", "KFFFL", "KYFIL", "AYFIL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
replicates = ["R1", "R2", "R3"]
base_dir = "data/contact_matrices"
output_dir = "results/heatmaps"
os.makedirs(output_dir, exist_ok=True)

amino_acid_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def get_custom_labels(sequence_name, original_labels):
    modified_labels = []
    if sequence_name == "KYYYL":
        y_count = 0
        for label in original_labels:
            if 'TYR' in label:
                y_count += 1
                modified_labels.append(f"Y{y_count}")
            else:
                modified_labels.append(amino_acid_map.get(label[:3].upper(), label))
    elif sequence_name == "KFFFL":
        f_count = 0
        for label in original_labels:
            if 'PHE' in label:
                f_count += 1
                modified_labels.append(f"F{f_count}")
            else:
                modified_labels.append(amino_acid_map.get(label[:3].upper(), label))
    else:
        modified_labels = [amino_acid_map.get(label[:3].upper(), label) for label in original_labels]
    return modified_labels

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)
    return LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )

# Plot style
plt.rcParams.update({'font.size': 10, 'font.family': 'Arial', 'axes.labelweight': 'bold'})

truncated_magma_cmap = truncate_colormap('magma', 0.1, 0.9)
heatmap_vmin, heatmap_vmax = 500, 2500

for seq in sequences:
    for rep in replicates:
        path = os.path.join(base_dir, seq, f"{rep}_contact_matrix.csv")
        if not os.path.exists(path):
            print(f"⚠️ Missing file: {path}")
            continue

        df = pd.read_csv(path, index_col=0)
        labels = get_custom_labels(seq, df.index.tolist())

        fig, ax = plt.subplots(figsize=(2.5, 2.5), dpi=300)
        sns.heatmap(df, ax=ax, cmap=truncated_magma_cmap,
                    cbar=False, vmin=heatmap_vmin, vmax=heatmap_vmax,
                    square=True, linewidths=0.5)

        ax.set_xticks(np.arange(len(labels)) + 0.5)
        ax.set_xticklabels(labels, fontsize=8, rotation=0)
        ax.set_yticks(np.arange(len(labels)) + 0.5)
        ax.set_yticklabels(labels, fontsize=8, rotation=0)
        ax.set_title(f"{seq} - {rep}", fontsize=10, fontweight='bold')

        cbar_ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
        sm = plt.cm.ScalarMappable(cmap=truncated_magma_cmap, norm=plt.Normalize(vmin=heatmap_vmin, vmax=heatmap_vmax))
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('Total Contacts', fontsize=10, fontweight="bold")

        plt.tight_layout(rect=[0, 0, 0.9, 1])
        output_path = os.path.join(output_dir, f"{seq}_{rep}_contact_heatmap.png")
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"✅ Saved: {output_path}")

print("🎨 Plotting complete. Check 'results/heatmaps/' for output PNGs.")
