import numpy as np
import seaborn as sns

# --- File Path and Naming Configuration ---
sequences = ["AYFIL", "KYFIL", "KAFIL", "KYAIL", "KYFIV", "KYFIA", "KYFAL", "KFFFL", "KYYYL"]
replicates = ["R1", "R2", "R3"]
# NOTE: Removed trailing slash from base_dir for cleaner string formatting
BASE_DIR = "/standard/BilodeauGroup/ShimantoRoy/Peptide/Peptide/KYFIL_replicated"

# --- System Name Generation ---
# The system name must combine the sequence and replicate (e.g., "AYFIL_R1")
# This creates a list of 27 unique system names (9 sequences * 3 replicates)
desired_order = [f"{seq}_{rep}" for seq in sequences for rep in replicates]


# --- MD Universe List Generation ---
# Generates the required list of tuples: (system_name, parm7_path, traj_path)
universes = []
for system_name in desired_order:
    # Example: system_name = "AYFIL_R1"
    
    # Path construction based on your file structure:
    # /standard/.../AYFIL_R1/parm7
    parm7_file = f"{BASE_DIR}/{system_name}/parm7"
    
    # /standard/.../AYFIL_R1/combined_trajectory.nc
    traj_file = f"{BASE_DIR}/{system_name}/combined_trajectory.nc"
    
    universes.append((system_name, parm7_file, traj_file))


# --- Other Configuration ---
# The color palette now accounts for the 27 total systems
color_palette = sns.color_palette("magma_r", len(desired_order))

# Trajectory frame parameters
start_frame = 10000
end_frame = 20001
step = 50

# RDF binning parameters
bin_width = 0.1
max_distance = 15.0
bins = np.arange(0, max_distance + bin_width, bin_width)