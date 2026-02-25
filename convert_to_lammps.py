# ==============================================================================
# SCRIPT: convert_to_lammps.py (DEBUG MODE)
# OBJECTIVE: Generate system.data with full error and progress tracing.
# ==============================================================================

import sys
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list

def log(msg):
    """Prints messages forcing buffer flush for instant viewing in the console."""
    print(msg, flush=True)

log("--- STARTING CONVERSION (DEBUG MODE) ---")

# 1. FILE LOADING
# ------------------------------------------------------------------------------
try:
    log("  [STEP 1] Loading structures...")
    atoms = read('system_solvated.pdb')
    
    # Box reference
    slab_ref = read('hematite_slab_final.pdb')
    ref_cell = slab_ref.get_cell()
    lx, ly = ref_cell[0, 0], ref_cell[1, 1]
    
    # Dynamic Lz
    z_coords = atoms.positions[:, 2]
    z_min, z_max = z_coords.min(), z_coords.max()
    lz = (z_max - z_min) + 5.0
    
    log(f"     -> System loaded: {len(atoms)} atoms.")
    log(f"     -> Box: {lx:.2f} x {ly:.2f} x {lz:.2f}")

except Exception as e:
    log(f"  [CRITICAL ERROR] Failed to load structures: {e}")
    sys.exit(1)

# 2. DEFINITIONS AND MEMORY ALLOCATION
# ------------------------------------------------------------------------------
MASSES = {1:55.845, 2:15.999, 3:15.999, 4:1.008, 5:15.999, 6:1.008, 
          7:22.990, 8:74.922, 9:15.999, 10:1.008}
CHARGES = {'fe':1.575, 'ob':-1.05, 'oh':-0.95, 'ho':0.425, 'ow':-0.8476, 
           'hw':0.4238, 'na':1.0, 'as':1.50, 'oa':-0.85, 'ha':0.45}

atom_types = np.zeros(len(atoms), dtype=int)
atom_charges = np.zeros(len(atoms))
molecule_ids = np.zeros(len(atoms), dtype=int)
bonds = []
angles = []

# Fixed Molecule IDs
MOL_SLAB = 1
MOL_ARSENATE = 2
MOL_SODIUM = 3
water_mol_counter = 4

# 3. PRE-CALCULATE NEIGHBORS (OPTIMIZED)
# ------------------------------------------------------------------------------
log("  [STEP 2] Calculating neighbor list (Optimization)...")
try:
    # Cutoff 1.25 detects O-H covalent bonds
    i_idx, j_idx = neighbor_list('ij', atoms, cutoff=1.25)
    
    # Convert to a list of lists for fast O(1) access
    adjacency_list = [[] for _ in range(len(atoms))]
    for k, atom_index in enumerate(i_idx):
        neighbor_index = j_idx[k]
        adjacency_list[atom_index].append(neighbor_index)
        
    log(f"     -> Connectivity map built for {len(atoms)} atoms.")

except Exception as e:
    log(f"  [ERROR] Failed in neighbor_list: {e}")
    sys.exit(1)

# 4. PASS 1: HEAVY ATOMS
# ------------------------------------------------------------------------------
log("  [STEP 3] Pass 1: Assigning Heavy Atoms...")
count_heavy = 0
try:
    for i, atom in enumerate(atoms):
        sym = atom.symbol
        molecule_ids[i] = MOL_SLAB # Default
        
        if sym == 'H': continue 

        # Use pre-calculated list
        neighbors = adjacency_list[i]
        neigh_syms = [atoms[n].symbol for n in neighbors]
        num_H = neigh_syms.count('H')

        if sym == 'Fe': 
            atom_types[i], atom_charges[i] = 1, CHARGES['fe']
        elif sym == 'Na': 
            atom_types[i], atom_charges[i] = 7, CHARGES['na']
            molecule_ids[i] = MOL_SODIUM
        elif sym == 'As': 
            atom_types[i], atom_charges[i] = 8, CHARGES['as']
            molecule_ids[i] = MOL_ARSENATE
        elif sym == 'O':
            # Arsenate detection (Distance-based backup)
            is_arsenate = False
            # Quick check only if it's not obvious
            if num_H != 2: 
                for n_idx, other in enumerate(atoms):
                    if other.symbol == 'As':
                        if np.linalg.norm(atom.position - other.position) < 2.0:
                            is_arsenate = True; break
            
            if is_arsenate:
                atom_types[i], atom_charges[i] = 9, CHARGES['oa']
                molecule_ids[i] = MOL_ARSENATE
            elif num_H == 2: # Water
                atom_types[i], atom_charges[i] = 5, CHARGES['ow']
                molecule_ids[i] = water_mol_counter
                water_mol_counter += 1
            elif num_H == 1: # Surface OH
                atom_types[i], atom_charges[i] = 3, CHARGES['oh']
            else: # Bulk O
                atom_types[i], atom_charges[i] = 2, CHARGES['ob']
        
        count_heavy += 1

    log(f"     -> {count_heavy} heavy atoms assigned.")

except Exception as e:
    log(f"  [ERROR] Failed in Pass 1 (Index {i}): {e}")
    sys.exit(1)

# 5. PASS 2: HYDROGENS
# ------------------------------------------------------------------------------
log("  [STEP 4] Pass 2: Assigning Hydrogens...")
count_h = 0
try:
    for i, atom in enumerate(atoms):
        if atom.symbol != 'H': continue

        neighbors = adjacency_list[i]
        if len(neighbors) == 0:
            log(f"     [WARN] Orphan hydrogen detected at ID {i+1}. Skipping.")
            continue

        parent_idx = neighbors[0]
        parent_type = atom_types[parent_idx]
        
        # Validate integrity
        if parent_type == 0:
            log(f"     [ERROR] H at {i+1} has parent {parent_idx+1} with no assigned type!")
            sys.exit(1)

        # Inherit Molecule ID
        molecule_ids[i] = molecule_ids[parent_idx]

        if parent_type == 5: # Water
            atom_types[i], atom_charges[i] = 6, CHARGES['hw']
        elif parent_type == 9: # Arsenate
            atom_types[i], atom_charges[i] = 10, CHARGES['ha']
        else: # Surface
            atom_types[i], atom_charges[i] = 4, CHARGES['ho']
            
        count_h += 1

    log(f"     -> {count_h} hydrogens assigned.")

except Exception as e:
    log(f"  [ERROR] Failed in Pass 2 (Index {i}): {e}")
    sys.exit(1)

# 6. PASS 3: TOPOLOGY (BONDS AND ANGLES)
# ------------------------------------------------------------------------------
log("  [STEP 5] Pass 3: Generating Bonds and Angles...")
try:
    for i in range(len(atoms)):
        # Progress trace every 500 atoms
        if i % 500 == 0: log(f"     ... processing topology for atom {i}/{len(atoms)}")

        my_type = atom_types[i]
        
        if my_type in [5, 3, 9]: # Ow, Oh, Oa
            h_neighbors = [n for n in adjacency_list[i] if atoms[n].symbol == 'H']
            
            # Water
            if my_type == 5:
                if len(h_neighbors) == 2:
                    h1, h2 = h_neighbors
                    bonds.append((1, i+1, h1+1))
                    bonds.append((1, i+1, h2+1))
                    angles.append((1, h1+1, i+1, h2+1))
                else:
                    log(f"     [WARN] Anomalous water at ID {i+1}: {len(h_neighbors)} H neighbors.")

            # Surface
            elif my_type == 3 and len(h_neighbors) == 1:
                bonds.append((2, i+1, h_neighbors[0]+1))
                
            # Arsenate
            elif my_type == 9 and len(h_neighbors) == 1:
                bonds.append((2, i+1, h_neighbors[0]+1))

    log(f"     -> Final topology: {len(bonds)} bonds, {len(angles)} angles.")

except Exception as e:
    log(f"  [ERROR] Failed in Pass 3 (Index {i}): {e}")
    sys.exit(1)

# 7. WRITING DATA FILE
# ------------------------------------------------------------------------------
log("  [STEP 6] Writing system.data...")
try:
    shift_vec = np.array([lx/2, ly/2, -z_min + 1.0]) 

    with open('system.data', 'w') as f:
        f.write("LAMMPS data file Generated by E2E Python Workflow\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n")
        f.write("0 dihedrals\n0 impropers\n\n")
        
        f.write("10 atom types\n")
        f.write("2 bond types\n")
        f.write("1 angle types\n\n")
        
        f.write(f"0.0 {lx:.6f} xlo xhi\n")
        f.write(f"0.0 {ly:.6f} ylo yhi\n")
        f.write(f"0.0 {lz:.6f} zlo zhi\n\n")
        
        f.write("Masses\n\n")
        for t in range(1, 11): f.write(f"{t} {MASSES[t]}\n")
        
        f.write("\nAtoms\n\n")
        # FULL FORMAT: ID Mol-ID Type Charge X Y Z
        for i in range(len(atoms)):
            p = atoms.positions[i] + shift_vec
            p[0] %= lx
            p[1] %= ly
            f.write(f"{i+1} {molecule_ids[i]} {atom_types[i]} {atom_charges[i]:.6f} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            
        f.write("\nBonds\n\n")
        for i, b in enumerate(bonds):
            f.write(f"{i+1} {b[0]} {b[1]} {b[2]}\n")
            
        f.write("\nAngles\n\n")
        for i, a in enumerate(angles):
            f.write(f"{i+1} {a[0]} {a[1]} {a[2]} {a[3]}\n")

    log("--- [SUCCESS] SCRIPT COMPLETED SUCCESSFULLY ---")

except Exception as e:
    log(f"  [CRITICAL ERROR] Failed while writing file: {e}")
    sys.exit(1)
