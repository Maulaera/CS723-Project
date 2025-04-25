import numpy as np
from Bio import PDB

def read_pdb(file):
    parser = PDB.PPBuilder()
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', file)
    return structure

def get_axis_points(structure):
    axis_points = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if 'CA' atom exists
                if 'CA' in residue:
                    axis_points.append(residue['CA'].get_coord())
                else:
                    # If 'CA' is missing, use the first available atom (e.g., 'N', 'C', 'O')
                    for atom in residue:
                        axis_points.append(atom.get_coord())
                        break  
    return np.array(axis_points)

def translate_points(points, translation_vector):
    return points + translation_vector

def rotate_points(points, rotation_matrix):
    return np.dot(points, rotation_matrix)

def align_axes(axis1_points, axis2_points):
    translation_vector = axis1_points[0] - axis2_points[0]

    translated_axis2 = translate_points(axis2_points, translation_vector)

    rotation_matrix = np.identity(3)  

    rotated_axis2 = rotate_points(translated_axis2, rotation_matrix)

    return rotated_axis2

def write_pdb(output_file, transformed_points):
    # Output the transformed PDB structure
    with open(output_file, 'w') as f:
        for i, point in enumerate(transformed_points):
            f.write(f"ATOM  {i+1:5d}  CA  ASN A {i+1:4d}    {point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}\n")

# Main function
def main():
    # Read input PDB files
    axis1_structure = read_pdb("helix1-axis.pdb")
    axis2_structure = read_pdb("helix2-axis.pdb")

    # Get axis points
    axis1_points = get_axis_points(axis1_structure)
    axis2_points = get_axis_points(axis2_structure)

    # Align the axes
    aligned_points = align_axes(axis1_points, axis2_points)

    # Write the aligned structure to a new PDB file
    write_pdb("aligned_helix2.pdb", aligned_points)

if __name__ == "__main__":
    main()
