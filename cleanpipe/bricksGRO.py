from cleanpipe import algelin
from cleanpipe import bricksStorage
from cleanpipe import bricksTOP
from cleanpipe import bricksFileSystem

import pandas as pd


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

            # Fixed-width fields (1-based columns in GROMACS spec):
            #  1–5:   residue number
            #  6–10:  residue name
            # 11–15:  atom name   (we ignore this column in your output)
            # 16–20:  atom number
            # 21–28:  x coordinate
            # 29–36:  y coordinate
            # 37–44:  z coordinate

            s_residue         = line[0:5].strip()
            s_structural_name = line[5:10].strip()
            s_atom_id         = line[15:20].strip()
            s_x               = line[20:28].strip()
            s_y               = line[28:36].strip()
            s_z               = line[36:44].strip()

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


def extract_all_dihedrals_from_gro(s_gro_file,s_top_file,s_molname):
    """
    will read the gro to find the real dihedrals as they are in the fro for the given molecule
    the top file will inform the list of dihedrals of the given molecule
    if there are several instantiations of that molecule. the dihedrals will come from the first instantiation


    example usage:

    s_gro_file = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.gro"
    s_top_file = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.top"
    cl.extract_all_dihedrals_from_gro(s_gro_file,s_top_file,'Protein_chain_A')
    """

    #parse the molecules in the top
    s_top_file_with_inclusions = bricksTOP.expand_includes_to_temp_file(s_top_file)
    dll_parsed_molecules = bricksTOP.parse_directives_inside_each_and_every_molecule(s_top_file_with_inclusions)
    bricksFileSystem.delete(s_top_file_with_inclusions)
    
    #get the dihedrals and store them in a df
    ll_dihedrals = dll_parsed_molecules[s_molname]['[ dihedrals ]']
    df_hihedrals = pd.DataFrame(ll_dihedrals)
    df_hihedrals = df_hihedrals.iloc[:, :4] # keep only the columns i,j,k,l, so to discart eventual parameters that might be set up manually
    df_hihedrals.columns=['i','j','k','l']
    
    
    #get the gobal id of the first atom in the molecule. this alows to calculate the global id for each atom of the moleule 
    dd_mols = bricksTOP.basic_infos_of_molecules(s_top_file)
    print(dd_mols)
    n_first_id = dd_mols[s_molname]['first_id']
    
    
    
    # Loop through dihedrals df row by row, to find the rows that define all phi and psi. then find their values "as is" according with the gro
    for index, row in df_hihedrals.iterrows():
    
        #get the local id and calculate the globa using the id of the first atom in the molecule
        global_i = int(row['i'])+(int(n_first_id)-1)
        global_j = int(row['j'])+(int(n_first_id)-1)
        global_k = int(row['k'])+(int(n_first_id)-1)
        global_l = int(row['l'])+(int(n_first_id)-1)
    
        #calculate de diheral angles
        df_hihedrals.at[index, 'diheral found in gro (1st molecule)'] = dihedral(s_gro_file, global_i, global_j, global_k, global_l)

    #this is an example of the df_hihedrals. the dihedrals are in the last column:
    
	#   i	j	k	l	functional	diheral found in gro (1st molecule)
    #0	2	1	5	6	9	        -61.084367
    #1	2	1	5	7	9	        58.161350
    #2	2	1	5	20	9	        180.000000
    #3	3	1	5	6	9	        58.752668

    
    return bricksStorage.df2ll(df_hihedrals) # return the information as a list of lists. be aware that the column names will be lost

