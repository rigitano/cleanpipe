import numpy as np

def find_new_atom_coord(atom1, atom2, atom3, distance, angle_in_plane_deg, angle_to_plane_deg):
    """
    
    code in python to find those x y and z coordinates using as inputs: 
    a list with the coordinates of atom 1, 
    a list with the coordinates of atom 2, 
    a list with the coordinates of atom 3, 
    a distance from atom 1, 
    an angle (referenced from the line formed from the the line between atoms 1 and 2) within the plane formed by the 3 atoms, 
    and a second angle (referenced from the line formed from the the line between atoms 1 and 2) within the plane formed by the 3 atoms

    """

    
    # Convert angles from degrees to radians
    angle_in_plane = np.radians(angle_in_plane_deg)
    angle_to_plane = np.radians(angle_to_plane_deg)
    
    # Vectors between atoms
    v12 = np.array(atom1) - np.array(atom2)
    v13 = np.array(atom3) - np.array(atom1)
    
    # Normal to the plane defined by the three atoms
    normal = np.cross(v12, v13)
    normal_unit = normal / np.linalg.norm(normal)

    # Unit vector along atom1 to atom2
    v12_unit = v12 / np.linalg.norm(v12)
    
    # Vector in the plane perpendicular to v12
    in_plane_perpendicular = np.cross(normal, v12)
    in_plane_perpendicular_unit = in_plane_perpendicular / np.linalg.norm(in_plane_perpendicular)
    
    # New atom position calculation
    direction_vector = (np.cos(angle_in_plane) * v12_unit + np.sin(angle_in_plane) * in_plane_perpendicular_unit)
                        
    direction_vector_unit = direction_vector / np.linalg.norm(direction_vector)
    
    # Final position calculation
    final_vector = np.cos(angle_to_plane) * direction_vector_unit + np.sin(angle_to_plane) * normal_unit
    final_position = np.array(atom1) + distance * final_vector

    return final_position