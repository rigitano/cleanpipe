import numpy as np
import math


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





def euclidian_distance(coord1: dict, coord2: dict) -> float:
    """
    Calculates the Euclidean distance between two 3D points.
    
    Args:
        coord1 (dict): A dictionary with keys 'x', 'y', 'z' representing the first coordinate.
        coord2 (dict): A dictionary with keys 'x', 'y', 'z' representing the second coordinate.
    
    Returns:
        float: The Euclidean distance between the two points.



    Example: 

    # Example usage:
    coord1 = {'x': 1, 'y': 2, 'z': 3}
    coord2 = {'x': 4, 'y': 6, 'z': 8}

    distance = distance(coord1, coord2)
    """

    return math.sqrt(
        (coord2['x'] - coord1['x']) ** 2 +
        (coord2['y'] - coord1['y']) ** 2 +
        (coord2['z'] - coord1['z']) ** 2
    )



def angle_in_middle_atom(coord1: dict, coord2: dict, coord3: dict) -> float:
    """
    Calculates the angle (in degrees) between two vectors defined by three points in 3D space.
    The second coordinate serves as the vertex of the angle.
    
    Args:
        coord1 (dict): The first coordinate.
        coord2 (dict): The vertex coordinate.
        coord3 (dict): The third coordinate.
    
    Returns:
        float: The angle in degrees between the two vectors.


    Example usage:
    coord1 = {'x': 1, 'y': 2, 'z': 3}
    coord2 = {'x': 4, 'y': 6, 'z': 8}
    coord3 = {'x': 7, 'y': 10, 'z': 12}

    angle = angle_between_vectors(coord1, coord2, coord3)


    """
    # Vector A (from coord2 to coord1)
    ax, ay, az = coord1['x'] - coord2['x'], coord1['y'] - coord2['y'], coord1['z'] - coord2['z']
    
    # Vector B (from coord2 to coord3)
    bx, by, bz = coord3['x'] - coord2['x'], coord3['y'] - coord2['y'], coord3['z'] - coord2['z']
    
    # Dot product of A and B
    dot_product = ax * bx + ay * by + az * bz
    
    # Magnitudes of A and B
    magnitude_a = math.sqrt(ax**2 + ay**2 + az**2)
    magnitude_b = math.sqrt(bx**2 + by**2 + bz**2)
    
    # Compute the angle in radians and convert to degrees
    if magnitude_a == 0 or magnitude_b == 0:
        return None  # Avoid division by zero
    
    cos_theta = dot_product / (magnitude_a * magnitude_b)
    cos_theta = max(-1, min(1, cos_theta))  # Clamp value to avoid floating point errors
    
    angle = math.degrees(math.acos(cos_theta))
    return angle




def dihedral_between_first3_and_last3(coord1: dict, coord2: dict, coord3: dict, coord4: dict) -> float:
    """
    Calculates the dihedral angle (in degrees) defined by four points in 3D space,
    following the GROMACS topology file convention.
    
    The angle lies between the second and third atom, with the first and fourth atoms as reference.
    
    Args:
        coord1 (dict): First reference point.
        coord2 (dict): First central point.
        coord3 (dict): Second central point.
        coord4 (dict): Second reference point.
    
    Returns:
        float: The dihedral angle in degrees.

    Example usage:
    coord1 = {'x': 1, 'y': 2, 'z': 3}
    coord2 = {'x': 4, 'y': 6, 'z': 8}
    coord3 = {'x': 7, 'y': 10, 'z': 12}
    coord4 = {'x': 10, 'y': 14, 'z': 16}

    dihedral = dihedral(coord1, coord2, coord3, coord4)

    """
    # Convert dictionaries to numpy arrays
    p1, p2, p3, p4 = np.array([coord1['x'], coord1['y'], coord1['z']]), \
                     np.array([coord2['x'], coord2['y'], coord2['z']]), \
                     np.array([coord3['x'], coord3['y'], coord3['z']]), \
                     np.array([coord4['x'], coord4['y'], coord4['z']])
    
    # Define bond vectors
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    
    # Normal vectors to planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    
    # Normalize vectors
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    b2 /= np.linalg.norm(b2)
    
    # Compute the dihedral angle
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, n2), b2)
    angle = np.degrees(np.arctan2(y, x))
    
    return angle


