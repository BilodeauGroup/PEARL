
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial',
    'axes.labelweight': 'bold',
    'axes.titlesize': 10,
    'axes.titleweight': 'bold',
    'legend.fontsize': 10,
    'legend.title_fontsize': 10
})

df = pd.read_csv("data/cluster_size_distribution.csv")
df = df[(df['Time_ns'] >= 200) & (df['Time_ns'] <= 400)]

def compute_weighted_avg(data):
    sizes = data['Cluster_Size'].values
    if len(sizes) > 0:
        weighted_avg = np.average(sizes, weights=sizes)
        error = sem(sizes)
    else:
        weighted_avg = 0
        error = 0
    return weighted_avg, error

desired_order = ["KYYYL", "KYFIL", "AYFIL", "KFFFL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
color_palette = sns.color_palette("magma", len(desired_order))

sequences = []
avg_cluster_sizes = []
sem_cluster_sizes = []

for seq in desired_order:
    reps = df[df['Sequence'] == seq]['Replicate'].unique()
    replicate_means = []
    for rep in reps:
        subset = df[(df['Sequence'] == seq) & (df['Replicate'] == rep)]
        avg, err = compute_weighted_avg(subset)
        if avg > 0:
            replicate_means.append(avg)
    if replicate_means:
        mean_across_reps = np.mean(replicate_means)
        propagated_sem = sem(replicate_means)
        sequences.append(seq)
        avg_cluster_sizes.append(mean_across_reps)
        sem_cluster_sizes.append(propagated_sem)

sequences = sequences[::-1]
avg_cluster_sizes = avg_cluster_sizes[::-1]
sem_cluster_sizes = sem_cluster_sizes[::-1]

plt.figure(figsize=(3, 3))
plt.barh(sequences, avg_cluster_sizes, xerr=sem_cluster_sizes, capsize=5,
         alpha=0.9, color=color_palette[:len(sequences)], edgecolor='None')
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_edgecolor('black')
plt.xlabel('Weighted Average Cluster Size', fontsize=10, weight='bold')
plt.ylabel('Peptide Sequence', fontsize=10, weight='bold')
plt.grid(False)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('cluster_Size_filtered.png', dpi=200, bbox_inches='tight')
plt.show()
