# PEARL – Peptide Ensemble Analysis & Research Library



PEARL is a modular analysis framework designed for studying self-assembling

pentapeptides using molecular dynamics (MD) simulations. This repository

contains computational workflows for quantifying peptide aggregation,

clustering, π–π interactions, residue contacts, and solvent structure. These

tools support structure–function insights into peptide hydrogel design,

sequence dependency, and emergent assembly behavior.



---



## 📁 Repository Structure



PEARL/

│── clustering_analysis/ # Cluster identification and time-series analysis

│── data/ # Trajectories, topology files, intermediate outputs

│── peptide_contact_analysis/ # Residue–specific contact calculations and plots

│── pi_pi_analysis_complete/ # Comprehensive π–π interaction detection pipeline

│── rdf_analysis/ # Radial distribution function workflows and outputs



yaml

Copy code



---



## 🔍 Module Descriptions



### **clustering_analysis/**

Contains scripts and notebooks for determining peptide cluster formation

across MD trajectories. Includes algorithms based on inter-residue distance

cutoffs, block averaging, and sequence comparisons.



Key outputs:

- Cluster count vs. time  

- Block-averaged cluster profiles  

- Sequence-dependent aggregation behavior  



---



### **data/**

Stores all necessary input and intermediate files for the analyses, including:

- Amber trajectory (`*.nc`) and topology (`parm7`) files  

- Preprocessed NumPy arrays  

- Contact matrices, block averages, RDF tables  



*Note:* Large files should ideally be excluded using `.gitignore`.



---



### **peptide_contact_analysis/**

Implements residue–residue contact detection based on a Lennard–Jones–derived

cutoff (4 Å, heavy atoms only). Generates time-averaged contact maps and

quantifies positional dependency of interactions across sequences.



Key features:

- Contact matrices  

- Heatmaps and publication-ready plots  

- Sequence comparison tools  



---



### **pi_pi_analysis_complete/**

Full pipeline for aromatic π–π interaction detection. Calculates the number and

identity of aromatic pairs (e.g., Phe–Tyr, Tyr–Tyr) that participate in stacking

interactions, including comparisons between OPC3 and OPC4 solvent models.



Outputs:

- Pair counts per frame  

- Time-averaged π–π interaction statistics  

- Bar plots with error bars  

- Per-sequence aromatic interaction profiles  



---



### **rdf_analysis/**

Tools for computing radial distribution functions for atom groups of interest,

supporting solvent structuring analysis around aromatic and hydrophobic motifs.



Outputs:

- g(r) curves  

- Integrated coordination numbers  

- Solvent structure comparison across sequences and solvent models  



---



## 🚀 Getting Started



### 1. Clone the repository

```bash

git clone <your-github-url>

cd PEARL



📊 Outputs & Analysis Goals

This repository supports analysis toward designing effective peptide
hydrogels by quantifying sequence-dependent self-assembly through:

Aggregation & clustering profiles

Heavy-atom contacts

Aromatic π–π interactions

Solvent structuring via RDF

Sequence comparisons in OPC3 vs OPC4 environments

These metrics help establish connections between peptide hydrophobicity,
aromaticity, positional context, and emergent assembly behavior.

👥 Authors & Acknowledgements

Shimanto Roy
Bilodeau Group & Lampe Group
University of Virginia

Special thanks to:

Dr. Camille Bilodeau

Dr. Kyle Lampe
for continuous support, feedback, and guidance.

📬 Contact

For questions or collaboration inquiries, please contact:

Shimanto Roy
Email: your email here
