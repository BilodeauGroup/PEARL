# RDF Analysis Project

This project performs radial distribution function (RDF) analysis on molecular dynamics trajectories using MDAnalysis.

## Structure
- `config.py`: Contains configuration parameters and system definitions.
- `main.py`: Computes RDFs and saves them as CSV files in the `data/` folder.
- `plot_rdf_data.py`: Reads RDF CSV files and generates plots.
- `data/`: Folder to store RDF CSV files.

## Usage
1. Run `main.py` to compute RDFs.
2. Run `plot_rdf_data.py` to visualize the RDFs.

## Requirements
Install dependencies using:
```bash
pip install -r requirements.txt
```
