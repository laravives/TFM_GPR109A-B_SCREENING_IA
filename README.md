# GPR109A/B Virtual Screening and AI Modeling

This repository contains the computational scripts developed for the TFM:

**"Comparative Computational Modeling of GPR109A/B Complexes with Small Molecules: Physics-Based Methods versus Artificial Intelligence Approaches”**

## Project Overview

The main objective of this project is to investigate the interaction of small molecules with the GPCR receptors GPR109A and GPR109B using different computational approaches.

The repository includes workflows related to:

- Molecular docking and virtual screening
- Molecular dynamics simulations and binding free energy calculations
- AI-based affinity prediction methods using Boltz-2 and FlowDock

The scripts were developed to facilitate reproducibility and workflow organization during the project.

---

## Repository Structure

```text
│
├── amber/
│   ├── prep.sh                 # System preparation
│   ├── min1.in                 # Energy minimization 1
│   ├── min2.in                 # Energy minimization 2
│   ├── heat.in                 # Heating stage
│   ├── eq.in                   # Equilibration
│   ├── prod.in                 # Production MD
│   ├── pmed.slurum             # SLURM submission script for running AMBER molecular dynamics workflows
│   └── mmpbsa.in               # MM/PBSA calculations
│
├── autodock/
│   ├── prepare_ligands.py
│   ├── prepare_receptors.py
│   └── run_screening.slurm
│
├── boltz/
│   ├── affinity.yaml
│   └── boltz_run.slurm
│
├── flowdock/
│   ├── prepare_csv.py
│   ├── prepare_receptor.sh
│   └── flowdock_run.slurm
│
├── .gitignore              
└── README.md

```

## Software and Computational Tools

Different software tools were employed throughout this project for structural preparation, virtual screening, affinity prediction, molecular dynamics simulations, and free energy analysis. The main tools used, following the approximate computational workflow, are listed below:

- **Python 3** for automation scripts, structural processing, and data analysis.
- **RDKit** for ligand manipulation and cheminformatics workflows.
- **Boltz-2** for protein–ligand affinity prediction using deep learning-based artificial intelligence models.
- **PDBFixer** for cleaning and preparing protein structures obtained from the Protein Data Bank.
- **Open Babel** for molecular file format conversion and additional ligand/receptor preparation.
- **Meeko** for preparing ligands and receptors in PDBQT format compatible with AutoDock workflows.
- **AutoDock-GPU** and **AutoGrid4** for GPU-accelerated molecular docking and affinity map generation.
- **Amber24** and **AmberTools24** for system preparation, molecular dynamics simulations, and MM/GBSA-MM/PBSA free energy calculations.
- **Antechamber** and **parmchk2** for ligand parameterization using the GAFF force field and AM1-BCC partial charge assignment.
- **tleap**, **cpptraj**, and **MMPBSA.py** for topology generation, trajectory analysis, and binding free energy calculations.
- **FlowDock** for protein–ligand complex generation and refinement using diffusion-based generative models.
- **DUDE-Z** for decoy compound generation.
