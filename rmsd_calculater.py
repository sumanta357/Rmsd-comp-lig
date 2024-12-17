import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import sys

def calculate_rmsd(reference_pdb, target_pdb, atom_selection="name CA"):
    try:
        # Load reference and target PDB files
        ref = mda.Universe(reference_pdb)
        target = mda.Universe(target_pdb)
        
        # Atom selection (debugging added)
        ref_atoms = ref.select_atoms(atom_selection)
        target_atoms = target.select_atoms(atom_selection)

        print(f"Reference atoms selected: {len(ref_atoms)}")
        print(f"Target atoms selected: {len(target_atoms)}")

        # Error handling: Check atom counts
        if len(ref_atoms) == 0 or len(target_atoms) == 0:
            print(f"Error: No atoms selected for {target_pdb}. Check atom names or selection syntax.")
            return None
        if len(ref_atoms) != len(target_atoms):
            print(f"Error: Atom counts mismatch between {reference_pdb} and {target_pdb}")
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
