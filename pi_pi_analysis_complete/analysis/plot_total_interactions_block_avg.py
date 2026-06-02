import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import pickle
import seaborn as sns
import numpy as np
import os
from matplotlib.patches import Patch

sequences = ["AYFIL", "KYFIL", "KAFIL", "KYAIL", "KYFIV", "KYFIA", "KYFAL", "KFFFL", "KYYYL"]
replicates = ["R1", "R2", "R3"]
block_size = 25

pickle_file = "data/pi_pi_interactions_data.pkl"
if not os.path.exists(pickle_file):
    print(f"Error: '{pickle_file}' not found.")
    exit()

try:
    with open(pickle_file, "rb") as f:
        pi_pi_data = pickle.load(f)
    print(f"✅ Loaded π-π interaction data from '{pickle_file}'")
except Exception as e:
    print(f"Error loading pickle file: {e}")
    exit()

desired_order = ["KYYYL", "KYFIL", "AYFIL", "KFFFL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
colors = sns.color_palette("magma_r", n_colors=len(desired_order))
color_map = dict(zip(desired_order, colors))

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

def block_average(data, block_size):
    if len(data) == 0:
        return np.array([])
    data_truncated = data[:(len(data) // block_size) * block_size]
    num_blocks = len(data_truncated) // block_size
    averaged_data = [np.mean(data_truncated[i * block_size:(i + 1) * block_size]) for i in range(num_blocks)]
    return np.array(averaged_data)

for seq in desired_order:
    plt.figure(figsize=(3, 3), dpi=300)
    data_plotted = False
    legend_handles = []

    for rep in replicates:
        label = f"{seq}_{rep}"
        if label in pi_pi_data:
            data = pi_pi_data[label]
            times_ns = np.array(data["times_ns"])
            parallel_counts = np.array(data["parallel_counts"])
            t_shaped_counts = np.array(data["t_shaped_counts"])
            total_counts = parallel_counts + t_shaped_counts
            seq_color = color_map[seq]

            truncated_len = (len(total_counts) // block_size) * block_size
            if truncated_len == 0:
                print(f"Skipping {label}: Not enough data.")
                continue

            times_ns_truncated = times_ns[:truncated_len]
            total_counts_truncated = total_counts[:truncated_len]

            block_min_counts = np.array([np.min(total_counts_truncated[i * block_size:(i + 1) * block_size]) for i in range(truncated_len // block_size)])
            block_max_counts = np.array([np.max(total_counts_truncated[i * block_size:(i + 1) * block_size]) for i in range(truncated_len // block_size)])
            block_mid_times_ns = np.array([np.mean(times_ns_truncated[i * block_size:(i + 1) * block_size]) for i in range(truncated_len // block_size)])

            plt.fill_between(block_mid_times_ns, block_min_counts, block_max_counts,
                             color=seq_color, alpha=0.1)

            averaged_total_counts = block_average(total_counts, block_size)
            averaged_times_ns = block_average(times_ns, block_size)

            if len(averaged_total_counts) > 0:
                line, = plt.plot(averaged_times_ns, averaged_total_counts,
                                 linestyle='-', linewidth=2,
                                 color=seq_color, alpha=0.9)
                legend_handles.append(line)
                data_plotted = True

    if data_plotted:
        plt.xlabel("Time (ns)", fontsize=10, fontweight="bold")
        plt.ylabel("Total π-π Interactions", fontsize=10, fontweight="bold")
        plt.title(seq, fontsize=10, fontweight="bold")
        original_range_patch = Patch(color=color_map[seq], alpha=0.1, label='Raw Data')
        full_legend_handles = [original_range_patch] + legend_handles
        plt.legend(handles=full_legend_handles, fontsize=6, loc='upper right')
        plt.tight_layout()
        plt.savefig(f"data/{seq}_total_pi_pi_interactions_averaged_shaded_raw.png")
        plt.close()
    else:
        print(f"No data for {seq}.")

print("✅ Time series plots saved.")