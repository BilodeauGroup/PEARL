"""
Compute contact matrices for peptide sequences using MDAnalysis.
Saves per-replicate CSV matrices into ./data/contact_matrices/{sequence}/.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import os

# Define sequence map
sequence_map = {
    "AYFIL": ["ALA", "TYR", "PHE", "ILE", "LEU"],
    "KYFIL": ["LYS", "TYR", "PHE", "ILE", "LEU"],
    "KAFIL": ["LYS", "ALA", "PHE", "ILE", "LEU"],
    "KYAIL": ["LYS", "TYR", "ALA", "ILE", "LEU"],
    "KYFIV": ["LYS", "TYR", "PHE", "ILE", "VAL"],
    "KYFIA": ["LYS", "TYR", "PHE", "ILE", "ALA"],
    "KYFAL": ["LYS", "TYR", "PHE", "ALA", "LEU"],
    "KFFFL": ["LYS", "PHE", "PHE", "PHE", "LEU"],
    "KYYYL": ["LYS", "TYR", "TYR", "TYR", "LEU"]
}

replicates = ["R1", "R2", "R3"]
base_dir = "/standard/BilodeauGroup/ShimantoRoy/Peptide/Peptide/KYFIL_replicated/"

# Constants
peptide_number = 5
cutoff = 4.0  # Angstroms
frame_start = 10001
frame_end = 20000
frame_step = 50

# Loop over all sequences and replicates
for seq_name, labels in sequence_map.items():
    for rep in replicates:
        print(f"Processing {seq_name} - {rep}")
        
        top_path = f"{base_dir}/{seq_name}_{rep}/parm7"
        traj_path = f"{base_dir}/{seq_name}_{rep}/combined_trajectory.nc"
        
        try:
            u = mda.Universe(top_path, traj_path, topology_format="PARM7")
        except Exception as e:
            print(f"Failed to load {seq_name} - {rep}: {e}")
            continue

        contact_prob = np.zeros((peptide_number, peptide_number))
        frame_count = 0

        for ts in u.trajectory[frame_start:frame_end+1:frame_step]:
            peptides = u.select_atoms("protein and not name H*")
            pept_fragments = peptides.fragments

            for i, frag_i in enumerate(pept_fragments):
                for j, frag_j in enumerate(pept_fragments):
                    if i == j:
                        continue
                    if np.min(distances.distance_array(frag_i.positions, frag_j.positions, box=u.dimensions)) < cutoff:
                        for s in range(peptide_number):
                            for l in range(peptide_number):
                                dist = distances.distance_array(
                                    frag_i.residues[s].atoms.positions,
                                    frag_j.residues[l].atoms.positions,
                                    box=u.dimensions
                                )
                                if np.min(dist) < cutoff:
                                    contact_prob[s, l] += 1
            frame_count += 1

        df = pd.DataFrame(contact_prob, index=labels, columns=labels)

        output_dir = f"data/contact_matrices/{seq_name}"
        os.makedirs(output_dir, exist_ok=True)
        df.to_csv(f"{output_dir}/{rep}_contact_matrix.csv")
        print(f"Saved: {output_dir}/{rep}_contact_matrix.csv")

print("✅ Contact matrix computation complete.")
