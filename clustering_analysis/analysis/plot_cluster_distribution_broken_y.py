
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'

df_sizes = pd.read_csv("data/cluster_size_distribution.csv")
df_sizes_filtered = df_sizes[(df_sizes["Time_ns"] >= 200) & (df_sizes["Time_ns"] <= 400)]
df_sizes_filtered["PeptideCount"] = df_sizes_filtered["Cluster_Size"]
replicate_peptide_counts = df_sizes_filtered.groupby(["Sequence", "Replicate", "Cluster_Size"])["PeptideCount"].sum().unstack(fill_value=0)

BREAK_LOWER = 1000
BREAK_UPPER = 2000
GLOBAL_MAX_Y_LIMIT = max(replicate_peptide_counts.values.max() * 1.1, BREAK_UPPER + 500)

desired_order = ["KYYYL", "KFFFL", "KYFIL", "AYFIL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
colors = sns.color_palette("magma_r", n_colors=len(desired_order))
color_map = dict(zip(desired_order, colors))

for seq in desired_order:
    if seq not in replicate_peptide_counts.index.get_level_values('Sequence').unique():
        print(f"No data found for sequence: {seq}. Skipping.")
        continue

    seq_replicates_data = replicate_peptide_counts.loc[seq]

    for replicate_num, data_series in seq_replicates_data.iterrows():
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, figsize=(2.5, 2.5), dpi=300,
                                                gridspec_kw={'height_ratios': [1, 3]})
        fig.subplots_adjust(hspace=0.05)

        x = data_series.index.astype(int)
        y = data_series.values
        current_color = color_map.get(seq, 'gray')

        ax_top.bar(x, y, color=current_color, edgecolor="none")
        ax_bottom.bar(x, y, color=current_color, edgecolor="none")

        ax_bottom.set_ylim(0, BREAK_LOWER)
        ax_top.set_ylim(BREAK_UPPER, GLOBAL_MAX_Y_LIMIT)

        ax_top.spines['bottom'].set_visible(False)
        ax_bottom.spines['top'].set_visible(False)
        ax_top.xaxis.tick_top()
        ax_top.tick_params(labeltop=False)
        ax_bottom.xaxis.tick_bottom()

        ax_top.set_xlim(0, 25)
        ax_bottom.set_xlim(0, 25)
        ax_bottom.set_xticks(range(0, 25, 3))
        ax_top.set_xticks(range(0, 25, 3))

        ax_top.set_title(f"{seq} -  {replicate_num}", fontsize=10, weight='bold')
        ax_bottom.set_xlabel("Cluster Size", fontsize=10, weight='bold')
        fig.supylabel("Weighted Count", fontsize=10, weight='bold', x=0.1)

        ax_top.tick_params(axis='y', labelsize=10)
        ax_bottom.tick_params(axis='y', labelsize=10)
        ax_bottom.tick_params(axis='x', labelsize=10)

        ax_top.grid(False)
        ax_bottom.grid(False)

        ax_top.axvline(x=5, color='red', linestyle='--', linewidth=1.5)
        ax_bottom.axvline(x=5, color='red', linestyle='--', linewidth=1.5)

        d = .015
        kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False)
        ax_top.plot((-d, +d), (-d, +d), **kwargs)
        ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax_bottom.transAxes)
        ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

        plt.tight_layout()
        file_name = f"cluster_distribution_broken_y_{seq}_Rep{replicate_num}.png"
        fig.savefig(file_name, bbox_inches='tight')
        plt.show(fig)
        print(f"✅ Plot saved as: {file_name}")
