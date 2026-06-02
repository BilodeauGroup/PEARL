import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

with open("data/pi_pi_interactions_data.pkl", "rb") as f:
    pi_pi_data = pickle.load(f)

desired_order = ["KYYYL", "KYFIL", "AYFIL", "KFFFL", "KAFIL", "KYAIL", "KYFAL", "KYFIV", "KYFIA"]
colors = sns.color_palette("magma_r", n_colors=len(desired_order))
peptide_colors = dict(zip(desired_order, colors))

start_time = 200.0
end_time = 400.0

grouped_totals = defaultdict(list)
grouped_replicate_counts = defaultdict(list)

for label, data in pi_pi_data.items():
    sequence = label.split('_')[0]
    times_ns = np.array(data["times_ns"])
    parallel_counts = np.array(data["parallel_counts"])
    t_shaped_counts = np.array(data["t_shaped_counts"])
    total_counts = parallel_counts + t_shaped_counts
    mask = (times_ns >= start_time) & (times_ns <= end_time)
    selected_counts = total_counts[mask]
    grouped_totals[sequence].append(np.sum(selected_counts))
    grouped_replicate_counts[sequence].extend(selected_counts)

results = {}
for seq in desired_order:
    replicate_sums = np.array(grouped_totals.get(seq, []))
    mean_total = np.mean(replicate_sums) if replicate_sums.size > 0 else 0
    sem = np.std(replicate_sums, ddof=1) / np.sqrt(len(replicate_sums)) if len(replicate_sums) > 1 else 0
    results[seq] = {"mean": mean_total, "sem": sem}

labels = list(results.keys())
totals = [results[seq]["mean"] for seq in labels]
sems = [results[seq]["sem"] for seq in labels]
bar_colors = [peptide_colors[seq] for seq in labels]

plt.figure(figsize=(3, 3), dpi=200)
plt.barh(labels[::-1], totals[::-1], xerr=sems[::-1], color=bar_colors[::-1], capsize=4)
plt.xlabel("Total π-π Interactions", fontsize=10, fontweight='bold')
plt.ylabel("Peptide Sequence", fontsize=10, fontweight='bold')
plt.tight_layout()
plt.savefig(f"data/pi_pi_total_interactions_fixed_sem_{int(start_time)}_{int(end_time)}.png", dpi=300)
plt.show()
