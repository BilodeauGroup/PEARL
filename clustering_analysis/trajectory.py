
import numpy as np
from MDAnalysis.analysis import distances
from .clustering import find_peptide_clusters

def process_trajectory(universe, peptide_selection, start_frame, end_frame, step):
    peptides = universe.select_atoms(peptide_selection)
    num_clusters = []
    cluster_data = []

    for frame in range(start_frame, end_frame, step):
        universe.trajectory[frame]
        time_ns = frame * 0.02

        binary_contact_matrix = np.zeros((len(peptides.fragments), len(peptides.fragments)), dtype=bool)

        for i, pep1 in enumerate(peptides.fragments):
            for j, pep2 in enumerate(peptides.fragments):
                if j > i:
                    pep1_heavy = pep1.select_atoms("not name H*")
                    pep2_heavy = pep2.select_atoms("not name H*")
                    dist_arr = distances.distance_array(pep1_heavy.positions, pep2_heavy.positions, box=universe.dimensions)

                    if np.min(dist_arr) < 4:
                        binary_contact_matrix[i, j] = binary_contact_matrix[j, i] = 1

        clusters = find_peptide_clusters(binary_contact_matrix)
        num_clusters.append(len(clusters))

        for cluster in clusters:
            cluster_data.append([time_ns, len(cluster)])

    return num_clusters, cluster_data
