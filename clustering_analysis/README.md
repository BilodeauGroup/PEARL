# Peptide Clustering Analysis

This project analyzes peptide clustering behavior using MDAnalysis and visualizes results using various analysis scripts.

## Structure
- `trajectory.py`: Processes trajectory and computes clusters.
- `clustering.py`: Cluster detection logic.
- `utils.py`: Frame selection utility.
- `config.py`: Configuration for sequences and replicates.
- `main.py`: Main script to run the analysis.
- `analysis/`: Contains plotting and statistical analysis scripts.

## Usage
Install dependencies:
```bash
pip install -r requirements.txt
```
Run the main analysis and generate CSV files:
```bash
python main.py
```
Run analysis scripts:
```bash
python analysis/plot_clusters_per_sequence.py
python analysis/plot_cluster_distribution_broken_y.py
python analysis/plot_cluster_area_with_sem.py
python analysis/plot_total_interactions_block_avg.py
python analysis/plot_total_interactions_sem_window.py
python analysis/plot_weighted_avg_cluster_size.py
```

## Output
- CSV files in `data/`
- Plots saved in working directory or `sequence_plots/`
