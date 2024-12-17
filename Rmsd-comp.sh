#!/bin/bash

# Function to print usage information
usage() {
    echo "Usage: $0 <reference.pdb> <molecule_folder> <output.csv>"
    echo " - <reference.pdb>: Path to the reference PDB file."
    echo " - <molecule_folder>: Folder containing target PDB files."
    echo " - <output.csv>: Output CSV file to save RMSD results."
    exit 1
}

# Check if the required inputs are provided
if [ "$#" -ne 3 ]; then
    usage
fi

# Input arguments
reference_pdb=$1
molecule_folder=$2
output_csv=$3

# Verify that inputs exist
if [ ! -f "$reference_pdb" ]; then
    echo "Error: Reference PDB file '$reference_pdb' not found."
    exit 1
fi

if [ ! -d "$molecule_folder" ]; then
    echo "Error: Molecule folder '$molecule_folder' not found."
    exit 1
fi

# Header for output CSV
echo "Molecule_Name,RMSD" > "$output_csv"

# Iterate through PDB files in the molecule folder
for target_pdb in "$molecule_folder"/*.pdb; do
    if [ ! -f "$target_pdb" ]; then
        echo "Warning: No PDB files found in $molecule_folder."
        continue
    fi
    
    # Extract molecule name
    molecule_name=$(basename "$target_pdb")

    # Calculate RMSD using MDAnalysis via an inline Python command
    rmsd=$(python3 -c "
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
try:
    ref = mda.Universe('$reference_pdb')
    target = mda.Universe('$target_pdb')
    ref_atoms = ref.select_atoms('all')
    target_atoms = target.select_atoms('all')
    if len(ref_atoms) != len(target_atoms):
        raise ValueError('Atom counts do not match between reference and target.')
    rmsd_analysis = RMSD(target_atoms, ref_atoms).run()
    rmsd_value = rmsd_analysis.rmsd[-1][2]  # Last RMSD value
    print(round(rmsd_value, 3))
except Exception as e:
    print('Error')
")

    # Handle errors and write to CSV
    if [[ "$rmsd" == "Error" ]]; then
        echo "$molecule_name,Error" >> "$output_csv"
        echo "Error calculating RMSD for $molecule_name"
    else
        echo "$molecule_name,$rmsd" >> "$output_csv"
        echo "RMSD for $molecule_name: $rmsd Å"
    fi
done

echo "RMSD calculations complete. Results saved in $output_csv."
