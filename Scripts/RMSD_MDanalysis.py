import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import pandas as pd

# Function to calculate RMSD
def calculate_rmsd(reference_file, trajectory_file, selection="backbone"):
    # Load the reference structure
    reference = mda.Universe(reference_file)

    # Load the trajectory
    trajectory = mda.Universe(reference_file, trajectory_file)

    # Select atoms for RMSD calculation (e.g., backbone atoms)
    selection_atoms = trajectory.select_atoms(selection)

    # Align trajectory to the reference structure
    align.AlignTraj(trajectory, reference, select=selection).run()

    # Calculate RMSD
    rmsd_analysis = rms.RMSD(
        trajectory, reference, select=selection, groupselections=[selection]
    )
    rmsd_analysis.run()

    # Store results in a DataFrame
    results_df = pd.DataFrame({
        "Time": rmsd_analysis.times,
        "RMSD": rmsd_analysis.rmsd[:, 2]  # Extracting the 3rd column (total RMSD)
    })

    # Save results to an Excel file
    excel_file = "rmsd_results.xlsx"
    results_df.to_excel(excel_file, index=False)
    print(f"Results saved to {excel_file}")

# Get user input
reference_structure = input("Enter the path to the reference structure (.gro): ")
md_trajectory = input("Enter the path to the MD trajectory (.trr): ")
selection = input("Enter the atom selection for RMSD calculation (e.g., backbone): ")

# Example usage
calculate_rmsd(reference_structure, md_trajectory, selection)
