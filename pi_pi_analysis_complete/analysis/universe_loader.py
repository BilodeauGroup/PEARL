import os
import MDAnalysis as mda

def load_universes(sequences, replicates, base_dir):
    universes = []
    for seq in sequences:
        for rep in replicates:
            parm7_file = os.path.join(base_dir, f"{seq}_{rep}", "parm7")
            traj_file = os.path.join(base_dir, f"{seq}_{rep}", "combined_trajectory.nc")
            if os.path.exists(parm7_file) and os.path.exists(traj_file):
                universes.append((f"{seq}_{rep}", parm7_file, traj_file))
    return universes
