
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial',
    'axes.labelweight': 'bold',
    'axes.titlesize': 10,
    'axes.titleweight': 'bold',
    'legend.fontsize': 6,
    'legend.title_fontsize': 10
})

df_clusters = pd.read_csv("data/num_clusters_over_time.csv")
desired_order = ["KYYYL", "KFFFL", "KYFIL", "AYFIL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
colors = sns.color_palette("magma_r", n_colors=len(desired_order))
markers = ['*', 's', '^']

def block_average(data, block_size):
    if len(data) < block_size:
        return np.mean(data)
    else:
        num_blocks = len(data) // block_size
        block_means = [np.mean(data[i * block_size: (i + 1) * block_size]) for i in range(num_blocks)]
        return np.array(block_means)

output_dir = "sequence_plots"
os.makedirs(output_dir, exist_ok=True)

for seq, color in zip(desired_order, colors):
    fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
    df_sequence = df_clusters[df_clusters["Sequence"] == seq]
    replicates = df_sequence["Replicate"].unique()

    for i, rep in enumerate(replicates):
        df_rep = df_sequence[df_sequence["Replicate"] == rep]
        time_ns = df_rep["Time_ns"] / 2
        num_clusters = df_rep["Num_Clusters"].values

        block_size = 30
        averaged_clusters = block_average(num_clusters, block_size)
        averaged_time_ns = block_average(time_ns, block_size)

        ax.plot(averaged_time_ns, averaged_clusters, label=f"Replicate {rep}", color=color, marker=markers[i], linestyle='-', markersize=4)

    ax.set_title(seq)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of Clusters")
    ax.set_xlim(0, 420)
    ax.set_ylim(0, 30)
    ax.legend()
    plt.tight_layout()
    filename = f"{seq}_clusters.png"
    filepath = os.path.join(output_dir, filename)
    fig.savefig(filepath)
    print(f"Saved {filepath}")
    plt.show()
