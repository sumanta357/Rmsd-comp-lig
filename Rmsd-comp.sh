#!/bin/bash

# Bash script to compare RMSD between a reference molecule and multiple molecules
# using Python and MDAnalysis for atom-name-based alignment.

# Input: Reference PDB file and folder containing target PDB files
REFERENCE_PDB="$1"
TARGET_DIR="$2"
OUTPUT_CSV="rmsd_results.csv"

# Check if inputs are valid
if [[ ! -f "$REFERENCE_PDB" ]]; then
    echo "Error: Reference PDB file not found!"
    exit 1
fi

if [[ ! -d "$TARGET_DIR" ]]; then
    echo "Error: Target directory not found!"
    exit 1
fi

# Header for the CSV output
echo "Molecule,RMSD (Å)" > "$OUTPUT_CSV"

# Python RMSD calculation function (embedded inline)
cat << EOF > rmsd_calculator.py
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import sys

def calculate_rmsd(reference_pdb, target_pdb, atom_selection="name CA"):
    try:
        # Load reference and target PDB files
        ref = mda.Universe(reference_pdb)
        target = mda.Universe(target_pdb)
        
        # Atom selection based on names (e.g., CA for alpha carbons)
        ref_atoms = ref.select_atoms(atom_selection)
        target_atoms = target.select_atoms(atom_selection)
        
        # Error handling: Check atom counts
        if len(ref_atoms) != len(target_atoms):
            print("Error: Atom counts do not match!")
            return None
        
        # Calculate RMSD
        rmsd = RMSD(target_atoms, ref_atoms).run()
        return round(rmsd.rmsd[-1][2], 3)
    except Exception as e:
        print(f"Error processing {target_pdb}: {e}")
        return None

if __name__ == "__main__":
    reference_pdb = sys.argv[1]
    target_pdb = sys.argv[2]
    rmsd_value = calculate_rmsd(reference_pdb, target_pdb)
    if rmsd_value is not None:
        print(rmsd_value)
    else:
        print("Error")
EOF

# Loop through all PDB files in the target directory
for pdb_file in "$TARGET_DIR"/*.pdb; do
    if [[ -f "$pdb_file" ]]; then
        # Get molecule name (filename without path)
        molecule_name=$(basename "$pdb_file")
        
        # Run Python RMSD calculator
        rmsd_value=$(python3 rmsd_calculator.py "$REFERENCE_PDB" "$pdb_file")
        
        # Check for errors
        if [[ "$rmsd_value" == "Error" ]]; then
            echo "$molecule_name,Error" >> "$OUTPUT_CSV"
        else
            echo "$molecule_name,$rmsd_value" >> "$OUTPUT_CSV"
        fi
    fi
done

# Cleanup
rm rmsd_calculator.py
echo "RMSD comparison complete! Results saved to $OUTPUT_CSV"
