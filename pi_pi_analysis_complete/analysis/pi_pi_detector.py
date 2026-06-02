import numpy as np
import MDAnalysis.lib.distances
from .ring_utils import find_all_rings, calculate_normal_vector

def find_pi_stacking_rings(ring_list, ag, distmax=5.0, distmin=2.0,
                            angle_threshold_parallel=30,
                            t_shaped_angle_lower=60,
                            t_shaped_angle_upper=120):
    ring_cogs = [ring.center_of_mass() for ring in ring_list]
    ring_cogs_array = np.array(ring_cogs)

    ix_rings = MDAnalysis.lib.distances.self_capped_distance(
        ring_cogs_array, max_cutoff=distmax, min_cutoff=distmin,
        box=ag.dimensions, return_distances=False
    )

    ring_couples_to_stack = [(a, b) for a, b in ix_rings
                             if ring_list[a].atoms.fragments != ring_list[b].atoms.fragments]

    ix_of_actually_pi_stacking = []
    ix_t_shaped_pi_stacking = []

    for a, b in ring_couples_to_stack:
        vector_a = calculate_normal_vector(ring_list[a])
        vector_b = calculate_normal_vector(ring_list[b])
        cos_val = np.dot(vector_a, vector_b) / (np.linalg.norm(vector_a) * np.linalg.norm(vector_b) + 1e-12)
        cos_val = np.clip(cos_val, -1.0, 1.0)
        alpha = np.rad2deg(np.arccos(cos_val))
        if alpha < angle_threshold_parallel:
            ix_of_actually_pi_stacking.append((a, b))
        elif t_shaped_angle_lower <= alpha <= t_shaped_angle_upper:
            ix_t_shaped_pi_stacking.append((a, b))

    return ix_of_actually_pi_stacking, ix_t_shaped_pi_stacking
