import MDAnalysis as mda
import pandas as pd
from clustering_analysis.trajectory import process_trajectory
from clustering_analysis.utils import determine_frame_selection
from clustering_analysis.config import sequences, replicates, base_dir

cluster_data = []
size_distribution = []

for seq in sequences:
    for rep in replicates:
        parm7_file = f"{base_dir}/{seq}_{rep}/parm7"
        traj_file = f"{base_dir}/{seq}_{rep}/combined_trajectory.nc"

        try:
            u = mda.Universe(parm7_file, traj_file, topology_format="PARM7")
            start_frame, end_frame = determine_frame_selection(seq, len(u.trajectory))

            num_clusters, cluster_sizes = process_trajectory(u, "protein", start_frame, end_frame, step=50)

            for i, num in enumerate(num_clusters):
                cluster_data.append([seq, rep, (start_frame + i * 50) * 0.02, num])

            for time_ns, size in cluster_sizes:
                size_distribution.append([seq, rep, time_ns, size])

            print(f"Processed: {seq}_{rep}")

        except Exception as e:
            print(f"Error processing {seq}_{rep}: {e}")

import os

# Create data folder if it doesn't exist
os.makedirs("data", exist_ok=True)

# Save CSV files into the data folder
pd.DataFrame(cluster_data, columns=["Sequence", "Replicate", "Time_ns", "Num_Clusters"]).to_csv("data/num_clusters_over_time.csv", index=False)
pd.DataFrame(size_distribution, columns=["Sequence", "Replicate", "Time_ns", "Cluster_Size"]).to_csv("data/cluster_size_distribution.csv", index=False)

print("Analysis complete. CSV files saved in 'data/' folder.")
