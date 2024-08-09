import numpy as np
import pandas as pd
import MDAnalysis as mda

def find_ions_within_distance(u, atom_numbers, ion_resname, distance_cutoff=4.0):
    user_atoms = u.select_atoms(f"bynum {' '.join(map(str, atom_numbers))}")
    ions = u.select_atoms(f"resname {ion_resname}")

    time_list = []
    ion_atoms_list = []

    for ts in u.trajectory:
        distances = mda.lib.distances.distance_array(user_atoms.positions, ions.positions, box=u.dimensions)
        ions_within_cutoff = np.sum(distances <= distance_cutoff, axis=1)

        time_list.append(ts.time)
        ion_atoms_within_cutoff = ions[np.where(distances <= distance_cutoff)[1]]
        ion_atoms_within_cutoff_numbers = ion_atoms_within_cutoff.indices
        ion_atoms_list.append(ion_atoms_within_cutoff_numbers)

    df_ions = pd.DataFrame(ion_atoms_list)
    df_ions.insert(0, "Time", time_list)
    df_ions.insert(1, "Atom_Count", np.sum(df_ions.iloc[:, 1:].notnull(), axis=1))

    return df_ions

def save_output(df, atom_number, output_file):
    df.to_excel(output_file, index=False)

topology_file = input("Enter path to the topology file (.dms): ")
trajectory_file = input("Enter path to the trajectory file (.trr): ")

u = mda.Universe(topology_file, trajectory_file)

atom_numbers = [int(x) for x in input("Enter atom numbers to analyze (comma-separated): ").split(",")]
ion_resname = input("Enter the ion residue name: ")

for atom_number in atom_numbers:
    df_ions = find_ions_within_distance(u, [atom_number], ion_resname)
    output_ions_file = f"ion_atoms_output_{atom_number}.xlsx"
    df_ions.columns = [f"Atom_{col}" if col != "Time" else col for col in df_ions.columns]
    save_output(df_ions, atom_number, output_ions_file)
