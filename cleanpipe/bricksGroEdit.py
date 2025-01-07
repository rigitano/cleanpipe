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

def coordinate_by_id(s_gro_file, s_atom_id, format='dictionary'):

    """
    Returns the coordinates of an atom with a certain id in a .gro file.
    the output can be in on othe the folowing 3 formats, as chosen in the format parameter:
    1-'list' a list (simple to understand), 
    2-'white_spaced_string_in_angstroms' a white-spaced string with coordinates converted to angstrons (to use in tcl scripts)
    3-'dictionary' a dictionary (intuitive to get values from).

    be careful, the id in the gro matches the top only for the first molecule on the gro, usually the protein.
    for the following molecules, a match of ids between the gro and top
    will require subtraction of the last id of the previous molecule.
    this matching is not considered by this function, the function will just return the coordinate for a given id

    """
    print("barrrrrr")
    print("s_atom_id")
    print(s_atom_id)

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
    