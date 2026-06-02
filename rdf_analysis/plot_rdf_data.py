import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from config import desired_order, color_palette

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial',
    'axes.labelweight': 'bold',
    'axes.titlesize': 10,
    'axes.titleweight': 'bold',
    'legend.fontsize': 6,
    'legend.title_fontsize': 10
})

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(7, 5), dpi=200, sharey=True)
axes = axes.flatten()

for i, system_name in enumerate(desired_order):
    file_path = f'data/{system_name}_rdf.csv'
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        columns = df.columns[1:]

        ax = axes[i]
        for j, col in enumerate(columns):
            style = ['-', '--', ':'][j % 3]
            ax.plot(df['Distance'], df[col], label=col, linewidth=1, linestyle=style, color=color_palette[i])
        ax.legend(fontsize=6, loc='upper right', frameon=False)
        ax.set_title(system_name, fontsize=10, weight='bold')
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_xlim(0, df['Distance'].max())
    else:
        axes[i].axis('off')

fig.text(0.5, 0.04, 'Distance (Å)', ha='center', va='center', fontsize=10, weight='bold')
fig.text(0.02, 0.5, 'Counts(r)', va='center', ha='center', rotation='vertical', fontsize=10, weight='bold')
plt.tight_layout(rect=[0.05, 0.08, 0.98, 0.95], h_pad=1.5, w_pad=1.5)
plt.savefig('cluster_rdf.png', dpi=200, bbox_inches='tight')
plt.show()
