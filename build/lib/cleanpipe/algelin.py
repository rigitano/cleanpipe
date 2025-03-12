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


def scale_coordinates(coord, scale_x, scale_y, scale_z):
    """
    Scales each coordinate by its respective factor.

    Parameters:
        coord (dict): A dictionary with keys 'x', 'y', 'z'.
        scale_x (float): Scaling factor for the x coordinate.
        scale_y (float): Scaling factor for the y coordinate.
        scale_z (float): Scaling factor for the z coordinate.

    Returns:
        dict: A new dictionary with the scaled coordinates.
    """
    return {
        'x': coord['x'] * scale_x,
        'y': coord['y'] * scale_y,
        'z': coord['z'] * scale_z
    }

def move_coordinates(coord, x, y, z):
    """
    Scales each coordinate by its respective factor.

    Parameters:
        coord (dict): A dictionary with keys 'x', 'y', 'z'.
        x (float): number to add to x coordinate.
        y (float): number to add to y coordinate.
        z (float): number to add to z coordinate.

    Returns:
        dict: A new dictionary with the scaled coordinates.
    """
    return {
        'x': coord['x'] + x,
        'y': coord['y'] + y,
        'z': coord['z'] + z
    }


def calc_intermediate_point(coord1, coord2, percent):
    """
    Calculate an intermediate point between two coordinates.
    
    Parameters:
      coord1 (dict): A dictionary with keys 'x', 'y', 'z' representing the first coordinate.
      coord2 (dict): A dictionary with keys 'x', 'y', 'z' representing the second coordinate.
      percent (float): A value between 0 and 100 indicating how far along the line the point should be.
                       0 returns coord1, 100 returns coord2, and 50 returns the midpoint.
    
    Returns:
      dict: A dictionary representing the intermediate point with keys 'x', 'y', 'z'.
      
    Raises:
      ValueError: If percent is not between 0 and 100 or if coord1 and coord2 don't have the same keys.
    """
    
    if not (0 <= percent <= 100):
        raise ValueError("percent must be between 0 and 100")
    
    # Ensure both dictionaries have the same keys
    if coord1.keys() != coord2.keys():
        raise ValueError("Both coordinates must have the same keys")
    
    fraction = percent / 100.0
    intermediate = {}
    
    for key in coord1:
        intermediate[key] = coord1[key] + (coord2[key] - coord1[key]) * fraction
        
    return intermediate