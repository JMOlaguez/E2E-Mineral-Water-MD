# E2E-Mineral-Water-MD: Automated Workflow for Hydrated Mineral Interfaces

![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![LAMMPS](https://img.shields.io/badge/LAMMPS-Compatible-blue)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)

An open-source, End-to-End (E2E) automated workflow to reliably simulate complex mineral-water interfaces (e.g., Hematite $\alpha$-Fe2O3) using the **CLAYFF force field** in **LAMMPS**.

## 🛑 The Problem This Solves
Simulating hydroxylated surfaces with flexible force fields often leads to severe numerical instabilities ("flying ice cube" effect), steric crashes, and topological errors (`Bond atoms missing` across periodic boundaries). This toolset automatically resolves the "hybrid topology" problem (distinguishing bulk vs. surface oxygens) and sanitizes coordinates for stable slab geometry calculations.

## ✨ Key Features
1. **Automated Topology Assignment:** A 3-step Python algorithm that automatically maps structural oxygens, surface hydroxyls, and water molecules without manual intervention.
2. **Coordinate Sanitization:** Unwraps coordinates and centers the slab to prevent explosive artifacts when applying `kspace_modify slab` (2D Ewald).
3. **Validated Stability:** Includes LAMMPS protocols utilizing a soft-restraint initial minimization and a stable Langevin thermostat integration (0.5 fs) that eliminates energetic drift.

## 📂 Repository Structure
* `/scripts`: Python codes for topology generation, coordinate sanitization, and workflow orchestration.
* `/lammps_inputs`: Production-ready LAMMPS scripts (`in.production`) parameterized for stability.
* `/benchmark_data`: Validated `system_final.data` coordinate files to test your local LAMMPS installation.

## 🚀 Quick Start
# Create a working directory
mkdir run_test && cd run_test
# Copy necessary files from the repository structure
cp ../structures/Hematita.cif .
cp ../lammps_inputs/in.* .
cp ../scripts/*.py .
# Run the orchestrator
python3 manuscript_e2e.py

To run the automated workflow on your local machine:

# 1. Clone the repository
git clone [https://github.com/JMOlaguez/E2E-Mineral-Water-MD.git](https://github.com/JMOlaguez/E2E-Mineral-Water-MD.git)
cd E2E-Mineral-Water-MD

# 2. Run the orchestrator
python3 scripts/E2E_Workflow.py
