# ==============================================================================
# PROJECT: Automated E2E Workflow for Mineral-Water Interfaces
# SCRIPT: build_system_E2E.py
# VERSION: 2.1 (Fix: ASE append/Atom syntax and full generation)
# DESCRIPTION: 
#   1. Loads Hematite and cleaves along the (001) plane.
#   2. "Prunes" exposed Fe atoms to guarantee an O-terminated surface.
#   3. Hydroxylates the surface (adds H to exposed O atoms).
#   4. Generates solvent molecules (H2O), contaminant (H2AsO4-), and counter-ion (Na+).
#   5. Writes the PACKMOL input file for solvation.
# ==============================================================================

import numpy as np
from ase import Atoms, Atom  # <--- Important: Atom (singular) for the loop
from ase.io import read, write
from ase.build import surface

# --- PHYSICOCHEMICAL PARAMETERS ---
DENSITY_WATER = 0.997   # g/cm3 (SPC/E @ 298K)
WATER_LAYER_Z = 40.0    # Height of the water layer (Angstroms)
BOND_OH = 0.97          # O-H bond length for surface hydroxyls
LAYER_TOL = 0.5         # Tolerance for identifying atomic layers along Z

print("\n--- STARTING CONSTRUCTION PROTOCOL ---")

# 1. Load Hematite
print("  [1] Loading Hematite.cif and cleaving (001) plane...")
bulk = read("Hematita.cif")
# Create a 3x3 supercell of the (001) surface, 4 layers thick, vacuum of 0 (we add it later)
slab = surface(bulk, (0, 0, 1), layers=4, vacuum=0.0)
slab = slab.repeat((3, 3, 1))

# 2. Prune Surface Fe
print("  [2] Identifying surface termination...")
# Get unique Z heights to identify atomic layers
z_positions = slab.positions[:, 2]

# Find maximum and minimum atomic heights
z_max = np.max(z_positions)
z_min = np.min(z_positions)

fe_to_delete = []
for i, atom in enumerate(slab):
    # If the highest/lowest atom is Fe, we mark it for deletion to ensure O-termination
    if atom.symbol == 'Fe':
        if abs(atom.position[2] - z_max) < LAYER_TOL or abs(atom.position[2] - z_min) < LAYER_TOL:
            fe_to_delete.append(i)

print(f"      - Pruned Fe atoms: {len(fe_to_delete)}")

# Delete in reverse order to preserve indices during the loop
for i in sorted(fe_to_delete, reverse=True):
    del slab[i]

# 3. Hydroxylate Surface
print("  [3] Hydroxylating exposed oxygens (O-H termination)...")
# Recalculate layers after pruning
z_positions = slab.positions[:, 2]
z_max = np.max(z_positions)
z_min = np.min(z_positions)

# Identify exposed O atoms (top and bottom)
exposed_o = []
for i, atom in enumerate(slab):
    if atom.symbol == 'O':
        if abs(atom.position[2] - z_max) < LAYER_TOL:
            exposed_o.append({'idx': i, 'pos': atom.position, 'type': 'top'})
        elif abs(atom.position[2] - z_min) < LAYER_TOL:
            exposed_o.append({'idx': i, 'pos': atom.position, 'type': 'bottom'})

print(f"      - Adding {len(exposed_o)} protons (H+)...")

# Add H atoms using the Atom class
for o in exposed_o:
    if o['type'] == 'top':
        h_pos = o['pos'] + np.array([0, 0, BOND_OH])
    else:
        h_pos = o['pos'] - np.array([0, 0, BOND_OH])
    
    new_h = Atom('H', position=h_pos)
    slab.append(new_h)

# Save the Final Slab
write("hematite_slab_final.pdb", slab)

print("  [4] Generating PDB topologies for molecules...")

# --- SOLVENT AND CONTAMINANTS ---
# Water (Approximate SPC/E geometry)
water = Atoms('OH2', positions=[[0, 0, 0], 
                                [0.757, 0.586, 0.0], 
                                [-0.757, 0.586, 0.0]])
write("water.pdb", water)

# Sodium Ion
na = Atoms('Na', positions=[[0,0,0]])
write("sodium.pdb", na)

# Arsenate Ion (H2AsO4-) - Approximate Td geometry
arsenate = Atoms('AsO4H2', positions=[
    [0.000, 0.000, 0.000],      # As
    [1.180, 1.180, 1.180],      # O1 (unprotonated)
    [-1.180, -1.180, 1.180],    # O2 (unprotonated)
    [1.180, -1.180, -1.180],    # O3 (protonated)
    [-1.180, 1.180, -1.180],    # O4 (protonated)
    [2.000, -2.000, -2.000],    # H1 (attached to O3)
    [-2.000, 2.000, -2.000]     # H2 (attached to O4)
])
write("arsenate.pdb", arsenate)

print("      - water.pdb, sodium.pdb, arsenate.pdb saved.")

# --- PACKMOL INPUT ---
print("  [5] Calculating volume for PACKMOL...")

# Calculate dimensions of the slab cell
cell = slab.get_cell()
lx, ly = cell[0,0], cell[1,1]
z_water_start = np.max(slab.positions[:, 2]) + 2.0
z_water_end = z_water_start + WATER_LAYER_Z

# Calculate number of water molecules based on target density
water_vol_cm3 = (lx * ly * WATER_LAYER_Z) * 1e-24 
mol_water = (water_vol_cm3 * DENSITY_WATER) / 18.015
num_water = int(mol_water * 6.022e23)

print(f"  [INFO] Simulation box: {lx:.2f} x {ly:.2f} Angstroms")
print(f"  [INFO] Water molecules to insert: {num_water}")

packmol_content = f"""
# Input generated automatically by build_system_E2E.py
tolerance 2.0
filetype pdb
output system_solvated.pdb

# 1. Hematite surface (fixed in space)
structure hematite_slab_final.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  centerofmass
end structure

# 2. Arsenate ion (H2AsO4-)
structure arsenate.pdb
  number 1
  # Restricted within the water box, avoiding boundaries
  inside box 2.0 2.0 {z_water_start + 4.0} {lx-2.0} {ly-2.0} {z_water_end - 4.0}
end structure

# 3. Sodium counter-ion (Na+)
structure sodium.pdb
  number 1
  inside box 2.0 2.0 {z_water_start + 4.0} {lx-2.0} {ly-2.0} {z_water_end - 4.0}
end structure

# 4. Solvent (SPC/E Water)
structure water.pdb
  number {num_water - 1}  # -1 because we inserted Arsenate/Na
  inside box 0.0 0.0 {z_water_start} {lx} {ly} {z_water_end}
end structure
"""

with open("packmol.inp", "w") as f:
    f.write(packmol_content)

print("  [6] File packmol.inp generated.")
print("--- STRUCTURAL CONSTRUCTION COMPLETED ---")
