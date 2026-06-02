
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import sem

df_sizes = pd.read_csv("data/cluster_size_distribution.csv")
df_sizes_filtered = df_sizes[(df_sizes["Time_ns"] >= 200) & (df_sizes["Time_ns"] <= 400)]

desired_order = ["KYYYL", "KFFFL", "KYFIL", "AYFIL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
reversed_order = desired_order[::-1]
colors = sns.color_palette("magma_r", n_colors=len(desired_order))
reversed_colors = colors[::-1]

replicate_areas = {}
for name in desired_order:
    replicate_areas[name] = {}
    for replicate in df_sizes_filtered["Replicate"].unique():
        filtered_sizes = df_sizes_filtered[
            (df_sizes_filtered["Sequence"] == name) &
            (df_sizes_filtered["Cluster_Size"] > 5) &
            (df_sizes_filtered["Replicate"] == replicate)
        ]["Cluster_Size"].tolist()
        replicate_areas[name][replicate] = sum(filtered_sizes)

mean_areas = []
sem_areas = []
for name in desired_order:
    areas = list(replicate_areas[name].values())
    if areas:
        mean_areas.append(np.mean(areas))
        sem_areas.append(sem(areas))
    else:
        mean_areas.append(0)
        sem_areas.append(0)

plot_df = pd.DataFrame({
    'Sequence': desired_order,
    'Mean_Area': mean_areas,
    'SEM_Area': sem_areas
})

plot_df = plot_df.set_index('Sequence').loc[reversed_order].reset_index()

plt.figure(figsize=(3, 3), dpi=200)
plt.barh(plot_df['Sequence'], plot_df['Mean_Area'], xerr=plot_df['SEM_Area'],
         color=reversed_colors, capsize=5, alpha=0.9, edgecolor='None')

plt.ylabel('Peptide Sequence', fontsize=10,fontweight="bold")
plt.xlabel('Area (Cluster Sizes > 5)', fontsize=10,fontweight="bold")
plt.xlim(0, max(plot_df['Mean_Area'] + plot_df['SEM_Area']) * 1.1)
plt.show()
