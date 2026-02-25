# ==============================================================================
# SCRIPT: sanitize_data.py
# OBJECTIVE: Center the system along the Z-axis using the MEDIAN to avoid 
#            Periodic Boundary Condition (PBC) wrapping errors, and add vacuum.
# ==============================================================================
import sys
import numpy as np  # Required for calculating the median

input_file = "system_heated.data"
output_file = "system_fixed.data"
PADDING_VACUUM = 80.0  # Amount of vacuum to add (Angstroms)

print(f"--- SANITIZING {input_file} (ROBUST MEDIAN METHOD) ---")

try:
    with open(input_file, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"ERROR: File {input_file} not found.")
    sys.exit(1)

# 1. Read original dimensions and find pointers
lz_orig = 0.0
header_end_line = 0
atoms_start_line = 0
box_line_idx = -1

for i, line in enumerate(lines):
    if "zlo zhi" in line:
        parts = line.split()
        zlo = float(parts[0])
        zhi = float(parts[1])
        lz_orig = zhi - zlo
        box_line_idx = i
        print(f"  > Original Z Box: {zlo:.2f} to {zhi:.2f} (Lz = {lz_orig:.2f})")
    if "Atoms" in line:
        header_end_line = i
        atoms_start_line = i + 2

# 2. New box: Add vacuum to Z
lz_new = lz_orig + PADDING_VACUUM
print(f"  > New Z Box (+ {PADDING_VACUUM} vacuum): 0.0 to {lz_new:.2f}")

# 3. READ ATOMS TO FIND THE CENTER OF MASS
z_unwrapped = []
atom_indices = []

current_idx = atoms_start_line
while current_idx < len(lines):
    line = lines[current_idx]
    parts = line.split()
    if not parts:
        current_idx += 1
        continue
    if parts[0].isalpha():
        break

    try:
        # Format: id mol type q x y z nx ny nz
        z = float(parts[6])
        nz = int(parts[9]) if len(parts) >= 10 else 0

        # Continuous real coordinate (unwrapped)
        z_real = z + (nz * lz_orig)

        z_unwrapped.append(z_real)
        atom_indices.append(current_idx)

    except ValueError:
        pass
    current_idx += 1

if not z_unwrapped:
    print("CRITICAL ERROR: No atoms were read.")
    sys.exit(1)

# 4. CALCULATE DISPLACEMENT (THE KEY CORRECTION)
# Instead of using min() or max(), we use median() to find the true center of the slab mass
z_center_mass = np.median(z_unwrapped)
box_center_target = lz_new / 2.0

# The shift determines how much we must move the slab so its center aligns with the box center
shift_z = box_center_target - z_center_mass

print(f"  > Center of Mass (Z): {z_center_mass:.2f}")
print(f"  > Target Box Center: {box_center_target:.2f}")
print(f"  > Applying Robust Shift: {shift_z:.2f}")

# 5. Rewrite Atoms
for i, line_idx in enumerate(atom_indices):
    line = lines[line_idx]
    parts = line.split()

    # Calculate new Z coordinate
    new_z = float(parts[6]) + shift_z
    
    # We strip out the image flags (nx ny nz) by only taking the first 7 columns
    # because the system is now unwrapped and centered in a larger box
    lines[line_idx] = f"{parts[0]} {parts[1]} {parts[2]} {parts[3]} {parts[4]} {parts[5]} {new_z:.5f}\n"

# 6. Write the rest of the file
lines[box_line_idx] = f"0.0 {lz_new:.5f} zlo zhi\n"

with open(output_file, 'w') as f:
    f.writelines(lines)

print(f"SUCCESS! {output_file} saved, centered, and padded with vacuum.")
