import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from rdf_analysis.config import universes, desired_order, bins, start_frame, end_frame, step

def filter_inter_fragment_distances(atomgroup1, atomgroup2, universe):
    dists = distances.distance_array(atomgroup1.positions, atomgroup2.positions, box=universe.dimensions)
    frag_ids_1 = atomgroup1.fragments
    frag_ids_2 = atomgroup2.fragments
    mask = np.array([[frag_ids_1[i] != frag_ids_2[j] for j in range(len(frag_ids_2))] for i in range(len(frag_ids_1))])
    return dists[mask]

os.makedirs("data", exist_ok=True)

for system_name in desired_order:
    universe_data = next((u for u in universes if u[0] == system_name), None)
    if universe_data:
        system_name, parm7_file, traj_file = universe_data
        universe = mda.Universe(parm7_file, traj_file, topology_format="PARM7")
        box_volume = np.prod(universe.dimensions[:3])

        peptides = universe.select_atoms("protein").residues
        peptide_sequences = {res.resid: res.resname for res in peptides}
        head_resname = peptide_sequences[min(peptide_sequences.keys())]
        tail_resname = peptide_sequences[max(peptide_sequences.keys())]

        head_atoms = universe.select_atoms(f'resname {head_resname} and name CA')
        tail_atoms = universe.select_atoms(f'resname {tail_resname} and name CA')

        head_head_rdf = np.zeros(len(bins) - 1)
        tail_tail_rdf = np.zeros(len(bins) - 1)
        head_tail_rdf = np.zeros(len(bins) - 1)
        num_frames = 0

        for ts in universe.trajectory[start_frame:end_frame:step]:
            num_frames += 1
            if len(head_atoms) > 1:
                dists = filter_inter_fragment_distances(head_atoms, head_atoms, universe)
                head_head_rdf += np.histogram(dists, bins=bins)[0]
            if len(tail_atoms) > 1:
                dists = filter_inter_fragment_distances(tail_atoms, tail_atoms, universe)
                tail_tail_rdf += np.histogram(dists, bins=bins)[0]
            if len(head_atoms) > 0 and len(tail_atoms) > 0:
                dists = filter_inter_fragment_distances(head_atoms, tail_atoms, universe)
                head_tail_rdf += np.histogram(dists, bins=bins)[0]

        shell_volumes = (4 / 3) * np.pi * (bins[1:]**3 - bins[:-1]**3)
        ideal_density = (len(head_atoms) + len(tail_atoms)) / box_volume

        head_head_rdf /= (shell_volumes * ideal_density * num_frames)
        tail_tail_rdf /= (shell_volumes * ideal_density * num_frames)
        head_tail_rdf /= (shell_volumes * ideal_density * num_frames)

        df = pd.DataFrame({
            'Distance': bins[:-1],
            f'{head_resname}-{head_resname}': head_head_rdf,
            f'{tail_resname}-{tail_resname}': tail_tail_rdf,
            f'{head_resname}-{tail_resname}': head_tail_rdf
        })

        df.to_csv(f'data/{system_name}_rdf.csv', index=False)
        print(f"Saved: data/{system_name}_rdf.csv")
