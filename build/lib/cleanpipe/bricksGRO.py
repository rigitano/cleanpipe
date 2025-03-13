from cleanpipe import algelin


def parse_gro(s_gro_file):
    """
    Parses a .gro file and extracts atom information.



    # how to access, for example, the x coordinate of the first item in the list
    list_of_dicts = cl.parse_gro('path_to_gro_file.gro')
    s_x = list_of_dicts[0]['x']


    ex:
    cl.parse_gro('path_to_gro_file.gro')

    """
    list_of_dicts = []
    with open(s_gro_file, 'r') as file:
        lines = file.readlines()
        if len(lines) < 3:
            raise ValueError("The .gro file must have at least three lines.")
        for line in lines[2:-1]:  # Skip the first two lines (header) and the last line (box size)

            line_parts = line.split()
            s_residue         = line_parts[0]
            s_structural_name = line_parts[1]
            s_atom_id         = line_parts[2]
            s_x               = line_parts[3]
            s_y               = line_parts[4]
            s_z               = line_parts[5]

            atom_info = {
                'residue': s_residue,
                'structural_name': s_structural_name,
                'id': s_atom_id,
                'x': s_x,
                'y': s_y,
                'z': s_z
            }

            list_of_dicts.append(atom_info)
                    
    return list_of_dicts

def coordinate_by_id(s_gro_file, n_atom_id, format='dictionary'):

    """
    Returns the coordinates of an atom with a certain id in a .gro file.
    the output can be in on othe the folowing 3 formats, as chosen in the format parameter:
    1-'list' a list (simple to understand), 
    2-'white_spaced_string_in_angstroms' a white-spaced string with coordinates converted to angstrons (to use in tcl scripts)
    3-'dictionary' a dictionary (intuitive to get values from).

    be careful, the id here is  global


    """

    s_atom_id = str(n_atom_id)
    list_of_dicts = parse_gro(s_gro_file)




    # Iterate through the list to find the atom with a certain id
    for atom_info in list_of_dicts:
        if atom_info['id'] == s_atom_id:
            x_coordinate = float(atom_info['x'])
            y_coordinate = float(atom_info['y'])
            z_coordinate = float(atom_info['z'])
            break   

    if format == 'list':
        return [x_coordinate, y_coordinate, z_coordinate]
    elif format == 'white_spaced_string_in_angstroms':
        return f"{(x_coordinate*10):.3f} {(y_coordinate*10):.3f} {(z_coordinate*10):.3f}"
    elif format == 'dictionary':
        return {'x': x_coordinate, 'y': y_coordinate, 'z': z_coordinate}
    


def distance(s_gro_file,n_atom1, n_atom2):
    """

    given the global ids of 2 atoms, returns the distance between them 

    example:
    cl.get_distance('npt.gro', 2, 5)

    """
    d_coord1 = coordinate_by_id(s_gro_file, n_atom1, format='dictionary')
    d_coord2 = coordinate_by_id(s_gro_file, n_atom2, format='dictionary')


    return algelin.euclidian_distance(d_coord1, d_coord2)


def angle(s_gro_file, n_atom1, n_atom2, n_atom3):
    """
    given the global ids of 3 atoms, returns the angle they define acording to gromacs convention

    example:
    cl.get_angle('npt.gro', 2, 5, 6)

    """
    d_coord1 = coordinate_by_id(s_gro_file, n_atom1, format='dictionary')
    d_coord2 = coordinate_by_id(s_gro_file, n_atom2, format='dictionary')
    d_coord3 = coordinate_by_id(s_gro_file, n_atom3, format='dictionary')


    return algelin.angle_in_middle_atom(d_coord1, d_coord2, d_coord3)


def dihedral(s_gro_file, n_atom1, n_atom2, n_atom3, n_atom4):
    """
    given the global ids of 4 atoms, returns the dihedral they define acording to gromacs convention
    the first 3 define he first plane. the last 3 define the second plane
    in practive, this define a torsion angle between the second and the third

    example:
    cl.get_dihedral('npt.gro', 2, 5, 6, 8)

    """
    d_coord1 = coordinate_by_id(s_gro_file, n_atom1, format='dictionary')
    d_coord2 = coordinate_by_id(s_gro_file, n_atom2, format='dictionary')
    d_coord3 = coordinate_by_id(s_gro_file, n_atom3, format='dictionary')
    d_coord4 = coordinate_by_id(s_gro_file, n_atom4, format='dictionary')


    return algelin.dihedral_between_first3_and_last3(d_coord1, d_coord2, d_coord3, d_coord4)


