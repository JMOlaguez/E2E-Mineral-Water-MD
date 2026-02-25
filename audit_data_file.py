import numpy as np
import sys

# --- CONFIGURATION ---
FILENAME = "system_final.data"  # Target data file
TOLERANCE_CHARGE = 1e-4         # Tolerance for charge neutrality
PADDING_Z = 0.5                 # Small error margin for box boundaries

def audit_lammps_data(filename):
    print(f"🔍 STARTING QUALITY AUDIT: {filename}")
    print("-" * 60)
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"❌ CRITICAL ERROR: File '{filename}' not found.")
        print("   Ensure this script and the .data file are in the same directory.")
        return

    # 1. READ HEADERS
    natoms = 0
    box_bounds = {'z': []}
    atom_section_idx = -1
    
    for i, line in enumerate(lines):
        clean = line.strip().lower()
        
        # Detect number of atoms
        if "atoms" in clean and len(line.split()) >= 2 and atom_section_idx == -1:
            try:
                # Match the header line (e.g., "3390 atoms")
                if line.split()[1] == "atoms":
                    natoms = int(line.split()[0])
            except: pass
            
        # Detect Z box bounds
        if "zlo zhi" in clean:
            box_bounds['z'] = [float(line.split()[0]), float(line.split()[1])]
            
        # Detect start of Atoms section
        if clean.startswith("atoms"):
            atom_section_idx = i
            print(f"   -> Atoms section detected at line {i+1}: '{line.strip()}'")
            break

    print(f"📋 Atoms declared in header: {natoms}")
    if 'z' in box_bounds and len(box_bounds['z']) == 2:
        print(f"📦 Simulation Box (Z): {box_bounds['z']}")
    print("-" * 60)

    # 2. READ ATOMS
    print("📖 Reading atom data...")
    atoms = []
    if atom_section_idx != -1:
        for i in range(atom_section_idx + 2, len(lines)):
            parts = lines[i].split()
            if not parts: 
                continue
            # End of atoms section (usually starts with Velos, Bonds, etc.)
            if parts[0].isalpha(): 
                break
                
            # Full format: atom-ID molecule-ID atom-type q x y z
            if len(parts) >= 7:
                try:
                    q = float(parts[3])
                    z = float(parts[6])
                    atoms.append({'q': q, 'z': z})
                except ValueError:
                    pass

    print(f"✅ Atoms read successfully: {len(atoms)} / {natoms}")
    if len(atoms) != natoms:
        print("❌ WARNING: The number of atoms read does not match the header.")
    print("-" * 60)

    # 3. AUDIT
    errors = 0
    
    # A. Charge Neutrality
    total_charge = sum(a['q'] for a in atoms)
    print(f"⚡ Total System Charge: {total_charge:.6f}")
    if abs(total_charge) > TOLERANCE_CHARGE:
        print("❌ FAILED: The system is NOT neutral. Danger for kspace PPPM.")
        errors += 1
    else:
        print("✅ PASSED: Charge neutrality confirmed (Neutral System).")

    # B. Box Integrity (Vacuum and Boundaries)
    if not atoms:
        print("❌ ERROR: Empty atom list.")
        return

    z_coords = [a['z'] for a in atoms]
    min_z, max_z = min(z_coords), max(z_coords)
    box_min, box_max = box_bounds['z'][0], box_bounds['z'][1]
    
    print(f"📏 Minimum Z coordinate: {min_z:.3f} (Box limit: {box_min})")
    print(f"📏 Maximum Z coordinate: {max_z:.3f} (Box limit: {box_max})")
    
    # Boundary leak check
    if min_z < (box_min - PADDING_Z) or max_z > (box_max + PADDING_Z):
        print("❌ FAILED: Atoms outside the simulation box.")
        errors += 1
    else:
        print("✅ PASSED: All atoms are within the vertical boundaries.")
        
    # Vacuum check (for Slab)
    # Verify top vacuum space
    vacuum_space = box_max - max_z
    print(f"🌌 Top vacuum space: {vacuum_space:.2f} Å")
    if vacuum_space > 10.0:
        print("✅ PASSED: Sufficient vacuum for kspace_modify slab.")
    else:
        print("❌ FAILED: Insufficient vacuum. LAMMPS might fail with periodic artifacts.")
        errors += 1
        
    print("-" * 60)
    if errors == 0:
        print(f"🚀 CONCLUSION: The file {FILENAME} is VALID for production runs.")
    else:
        print(f"🛑 CONCLUSION: The file {FILENAME} has ERRORS. Do not simulate!")

if __name__ == "__main__":
    audit_lammps_data(FILENAME)
