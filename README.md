# MD Contact & Interaction Energy Analysis Pipeline 🔬

This repository contains a Python-based workflow designed to analyze molecular dynamics (MD) trajectories. The pipeline identifies key interacting residues between two molecular selections, quantifies their contact frequency, and sets up calculations to determine their non-bonded interaction energies (Electrostatics and van der Waals) using cpptraj.

The primary goal is to automate the process of pinpointing residues that maintain significant contact throughout a simulation and then calculating their specific energetic contributions to the interaction.

## Features

- Automated Workflow: A single command executes a multi-step analysis pipeline.
- Contact Analysis: Uses getcontacts[https://github.com/getcontacts/getcontacts] to identify and quantify residue-residue contacts.
- Visual Mapping: Maps contact frequencies onto the B-factor column of a PDB file for easy visualization.
- Intelligent Filtering: Identifies residues with the highest contact frequencies.
- Energy Calculation: Generates and optionally runs a cpptraj script to calculate LIE.
- Data Summarization: Parses cpptraj output into summary CSV tables.
- Configuration-Driven: All parameters are controlled via a JSON config file.

## Prerequisites

- Python 3
- AmberTools (cpptraj and ambpdb)
- Python libraries: pandas (install via `pip install pandas`)
- getcontacts[https://github.com/getcontacts/getcontacts]

## 🚀 Usage

1. Prepare a `config.json` file with paths and parameters.
2. Run the pipeline with:

    python full_pipeline.py config.json

## 📜 Workflow Summary

1. Calculate dynamic contacts
2. Compute contact frequencies
3. Extract reference PDB with ambpdb
4. Map contact frequencies into B-factor field
5. Identify residues with high contact frequency
6. Generate cpptraj script for LIE calculation
7. Run cpptraj and summarize interaction energies (if enabled)

## 📂 Output Files

- `dynamic_contacts.tsv`
- `contact_frequencies.tsv`
- `reference.pdb`
- `*_bfac.pdb`
- `/energies/lie_script.in`
- `/energies/int_*.dat`
- `/energies/ele_summary.csv`
- `/energies/vdw_summary.csv`
"""
