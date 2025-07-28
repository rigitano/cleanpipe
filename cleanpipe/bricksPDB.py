import PeptideBuilder
from Bio.PDB import PDBParser, PDBIO, Atom, Residue, Chain, Model, Structure
import numpy as np
import Geometry
from cleanpipe import bricksFileSystem
from cleanpipe import bricksPeptide
import subprocess
import string



def download_and_clean_pdb(s_molecule_name):
    """
    usage example:
    cl.download_and_clean_pdb("1aki")

    
    """

    #get pdb from portal
    subprocess.run(f"wget https://files.rcsb.org/download/{s_molecule_name}.pdb" , shell=True, check=True)

    #remove water
    subprocess.run(f"grep -v 'HOH' {s_molecule_name}.pdb > {s_molecule_name}_temp.pdb" , shell=True, check=True)
    bricksFileSystem.delete(f"{s_molecule_name}.pdb")
    subprocess.run(f"mv {s_molecule_name}_temp.pdb {s_molecule_name}.pdb" , shell=True, check=True)


def create_peptide(s_outName, s_aminoacids, l_phi, l_psi_im1, s_nTerminusCAP, s_cTerminusCAP ):

    """
    
    example:
    cl.create_peptide("pepticat2.pdb","AAAAAA", [-57.8,-57.8,-57.8,-57.8,-57.8,-57.8], [-47.0,-47.0,-47.0,-47.0,-47.0,-47.0], "ACE","NME" )
    or
    cl.create_peptide("pepticat2.pdb","AAAAAA", [-57.8,-57.8,-57.8,-57.8,-57.8,-57.8], [-47.0,-47.0,-47.0,-47.0,-47.0,-47.0], "","" )
    """

     #################################### create aminoacid chain ###################################

    # Add the rest of the amino acids to the peptide
    for i in range(0,len(s_aminoacids)):

        #get current aminoacid letter code, phi and psi_im1
        current_aminoacid = s_aminoacids[i]
        current_phi       = l_phi[i]
        current_psi_im1   = l_psi_im1[i]

        #define the geometry of the current aminoacid
        current_aa_geometry = PeptideBuilder.Geometry.geometry(current_aminoacid)
        current_aa_geometry.phi = current_phi
        current_aa_geometry.psi_im1 = current_psi_im1 #angle of the immediately preceding residue (where "im1" denotes "i minus 1")

        #insert the current aminoacid in peptide
        if i == 0:
            peptide = PeptideBuilder.initialize_res(current_aa_geometry) # Initialize the peptide with the first amino acid
        else:
            PeptideBuilder.add_residue(peptide, current_aa_geometry) # Add current aminoacid to previouly constructed peptide


    #################################### add termini ###################################
    if s_nTerminusCAP == "ACE":
        bricksPeptide.add_ACE_to_Nterminus(peptide)

    if s_cTerminusCAP == "NME":
        bricksPeptide.add_NME_to_Cterminus(peptide)


    #################################### create system. (ps this will add hydrogens) ###################################

    # Save temporary pdb file of the peptide
    io = PDBIO()
    io.set_structure(peptide)
    io.save(s_outName)






def insert_new_molecule_into_pdb(s_input_pdb_file, s_output_pdb_file, d_new_molecule_atoms):
    """
    Insert a new molecule into a pre-existing PDB file as a new chain.

    Parameters:
        input_pdb_file (str): Path to the input PDB file.
        output_pdb_file (str): Path to save the modified PDB file.
        new_molecule_atoms (list of dict): Each dict contains atom infos:

    # that molecule is something like this: 
    #[{'structural_name': 'XA','element_name': 'X', 'coord': [-3.346065214951231, 0.89657547216805345, 2.449489742783178], 'residue_name': 'TRS', 'residue_id': 1}, 
    # {'structural_name': 'XB','element_name': 'X', 'coord': [0.8965754721680534, -3.3460652149512313, 2.449489742783178], 'residue_name': 'TRS', 'residue_id': 1}, 
    # {'structural_name': 'XC','element_name': 'X', 'coord': [0.2347744239847329, -5.2938793894289373, 4.239847893938439], 'residue_name': 'TRS', 'residue_id': 1}, 
    # {'structural_name': 'XD','element_name': 'X', 'coord': [3.3460652149512313, -0.8965754721680534, -2.44948974278317], 'residue_name': 'TRS', 'residue_id': 1}]


    cl.insert_new_molecule_into_pdb("example.pdb", "modified_example.pdb", new_molecule_atoms) 
    """

    print("IN FUNCTION insert_new_molecule_into_pdb")
    #print(d_new_molecule_atoms)

    # Parse the existing PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', s_input_pdb_file)

    # Get the first model or create a new one if needed
    model = structure[0] if len(structure) > 0 else Model.Model(0)

    # Determine existing chain IDs to avoid collisions
    existing_chain_ids = {chain.id for chain in model.get_chains()}

    # Determine the chain ID for the new molecule
    new_chain_id = next((char for char in string.ascii_uppercase if char not in existing_chain_ids), None)
    if new_chain_id is None:
        raise ValueError("No unique chain ID available.")

    # Add a new chain for the molecule
    new_chain = Chain.Chain(new_chain_id)
    model.add(new_chain)

    # Determine the maximum residue ID in the whole structure to avoid collisions
    max_serial_number = max((atom.serial_number for atom in model.get_atoms()), default=0)

    # Insert residues and atoms into the new chain
    current_serial_number = max_serial_number + 1


    for current_atom in d_new_molecule_atoms:

        #get residue name and ensure that it is at most 3 characters to fit PDB conventions
        residue_name = current_atom['residue_name']
        residue_name = residue_name.ljust(3)[:3] 

        # Determine if the residue is a heteroatom based on its name
        isHeteroatom = residue_name not in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

        #check if the residue is already defined before
        #if not, add it to the chain
        #if yes, get the residue id
        # Check if the residue is already defined
        existing_residue = None
        for residue in new_chain.get_residues():
            if residue.id[1] == current_atom['residue_id']:
                #print(f"Residue {residue.id} already exists")
                existing_residue = residue
                break


        if existing_residue is None:
            # Create and add a new residue

            #define residue info (het_flag, resseq, icode)
            if isHeteroatom:
                residue_id = ('H_', current_atom['residue_id'], ' ')  
            else:
                residue_id = (' ', current_atom['residue_id'], ' ')

            #create and add residue
            residue = Residue.Residue(residue_id, residue_name, '')
            new_chain.add(residue)
        else:
            residue = existing_residue


        # Create new atom
        atom = Atom.Atom(
            current_atom['structural_name'].upper().ljust(4),  # Atom name (e.g., 'C1', 'N1', 'O1') # Ensure atom name is 4 characters (padded or truncated)
            np.array(current_atom['coord'], dtype=float),  # Coordinates as a numpy array
            1.0,  # B-factor (optional, default 1.0)
            1.0,  # Occupancy (optional, default 1.0)
            ' ',  # Alternate location indicator
            str(current_atom['structural_name']).strip(),  # Full atom name
            int(current_serial_number),  # Automatically assigned serial number
            current_atom['element_name'].upper() # the element name (e.g., 'C', 'N', 'O')
        )

        #print(f"Atom: {atom.get_name()}, coord: {atom.get_coord()}, Residue: {residue_name}, Residueid: {current_atom['residue_id']}, Serial Number: {atom.serial_number}")

        # Add atom to the residue
        residue.add(atom)

        # Increment residue ID for the next residue
        current_serial_number += 1
        # Increment residue ID for the next residue

    # Write the modified structure to the output file
    io = PDBIO()
    io.set_structure(structure)
    io.save(s_output_pdb_file)
    print(f"Modified PDB saved to {s_output_pdb_file}")



def add_truss(s_pdb_file, s_out_pdb_file, p1, p2, n_square_size = 3):
    """

    cl.add_truss("pepticat6_with_atom.pdb","pepticat6_with_atom_and_truss.pdb", [0,0,0], [10,10,10])

    """

    print("IN FUNCTION add_truss")
    

    # Calculate the direction vector and the total distance
    direction = np.array(p2, dtype=float) - np.array(p1, dtype=float)
    total_distance = np.linalg.norm(direction)
    
    # Normalize the direction vector
    direction = direction / total_distance
    
    # Calculate the number of points based on the square size
    n = int(total_distance // n_square_size)
    
    # Generate the points separated by the square size
    points = [np.array(p1, dtype=float) + i * n_square_size * direction for i in range(n + 1)]
    
    # Add the last point (p2)
    #points.append(np.array(p2, dtype=float))
    
    # Calculate two perpendicular vectors to the direction
    if direction[0] != 0 or direction[1] != 0:
        perp_vector1 = np.cross(direction, [0, 0, 1])
    else:
        perp_vector1 = np.cross(direction, [0, 1, 0])
    perp_vector1 = perp_vector1 / np.linalg.norm(perp_vector1)
    perp_vector2 = np.cross(direction, perp_vector1)
    
    # Generate the squares
    squares = []
    for point in points:
        # Create the square's vertices centered at the point
        square_vertices = []
        for dx, dy in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
            vertex = point + (dx * n_square_size / 2) * perp_vector1 + (dy * n_square_size / 2) * perp_vector2
            square_vertices.append(vertex)
        squares.append(square_vertices)
    
    # Collect vertices and edges
    vertices = []
    edges = []
    
    # Flatten the list of squares and generate vertices
    for square in squares:
        for vertex in square:
            vertices.append(vertex.tolist())
    
    # Generate edges within each square and between consecutive squares
    num_vertices_per_square = 4
    for i in range(len(squares)):
        # Edges within the square
        for j in range(num_vertices_per_square):
            edges.append([i * num_vertices_per_square + j, i * num_vertices_per_square + (j + 1) % num_vertices_per_square])
        
        # Edges between consecutive squares
        if i < len(squares) - 1:
            for j in range(num_vertices_per_square):
                for k in range(num_vertices_per_square):
                    edges.append([i * num_vertices_per_square + j, (i + 1) * num_vertices_per_square + k])
    
    # Calculate distances for each edge
    distances = []
    for edge in edges:
        v1 = np.array(vertices[edge[0]])
        v2 = np.array(vertices[edge[1]])
        distance = np.linalg.norm(v1 - v2)
        truncated_distance = int(distance * 1000) / 1000.0 #make it 3 decimal
        distances.append((edge[0], edge[1], truncated_distance))

    #inset the coordinates in a pdb file. will create a representation of the molecule that is list of dictionaries
    # than that I can pass that representation to function that inserts the molecule in the pdb
    #here is an example of such a list
    #[{'structural_name': 'XA','element_name': 'X', 'coord': [-3.346065214951231, 0.89657547216805345, 2.449489742783178], 'residue_name': 'TRS', 'residue_id': 1}, 
    # {'structural_name': 'XB','element_name': 'X', 'coord': [0.8965754721680534, -3.3460652149512313, 2.449489742783178], 'residue_name': 'TRS', 'residue_id': 1}, 
    # {'structural_name': 'XC','element_name': 'X', 'coord': [0.2347744239847329, -5.2938793894289373, 4.239847893938439], 'residue_name': 'TRS', 'residue_id': 1}, 
    # {'structural_name': 'XD','element_name': 'X', 'coord': [3.3460652149512313, -0.8965754721680534, -2.44948974278317], 'residue_name': 'TRS', 'residue_id': 1}]


    l_new_molecule_atoms = []
    #this loop will insert the four vertices of the square. each square is considered a residue. the for loop is repetead until all the squares are added
    residue_count=1
    for i in range(0, len(vertices), 4):
        l_new_molecule_atoms.append({'structural_name': f"XA",'element_name': 'X', 'coord': vertices[i+0], 'residue_name': 'TRS', 'residue_id': residue_count})
        l_new_molecule_atoms.append({'structural_name': f"XB",'element_name': 'X', 'coord': vertices[i+1], 'residue_name': 'TRS', 'residue_id': residue_count})
        l_new_molecule_atoms.append({'structural_name': f"XC",'element_name': 'X', 'coord': vertices[i+2], 'residue_name': 'TRS', 'residue_id': residue_count})
        l_new_molecule_atoms.append({'structural_name': f"XD",'element_name': 'X', 'coord': vertices[i+3], 'residue_name': 'TRS', 'residue_id': residue_count})
        #print(vertices[i])
        residue_count+=1

    #at last, we insert the new molecule in the pdb
    insert_new_molecule_into_pdb(s_pdb_file, s_out_pdb_file, l_new_molecule_atoms)


def insert_residue_into_chain(s_input_pdb_file,s_output_pdb_file,s_chain,n_residue,s_new_residue_name, d_new_residue_atoms,s_replace_or_displace='displace'):
    """


    example:
    d_hi =[{'structural_name': 'XI','element_name': 'X', 'coord': [1, 2, 3]}, 
        {'structural_name': 'XJ','element_name': 'X', 'coord': [4, 5, 6]}, 
        {'structural_name': 'XK','element_name': 'X', 'coord': [7, 8, 9]}, 
        {'structural_name': 'XL','element_name': 'X', 'coord': [10, 11, 12]}
      ]
    cl.insert_residue_into_chain("pepticat6_with_atom_and_truss.pdb","pepticat6_with_atom_and_truss_modified.pdb",'A',3,"BOB", d_hi,s_replace_or_displace='displace')
    """

    print("IN FUNCTION insert_residue_into_chain")

    #define the protein/heteroatom flag. this will be used latter. biopython will require this weird flag
    if s_new_residue_name not in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']:
        het_flag = 'H_'
    else:
        het_flag = ' '


    # Parse the existing PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', s_input_pdb_file)

    # Get the first model or create a new one if needed
    model = structure[0] if len(structure) > 0 else print('CLEAN PIPE MESSAGE : ATTENTION there are several models in the pdb file, it is the first that will be edited')

    # Determine the maximum residue ID in the whole structure to avoid collisions
    #max_serial_number = max((atom.serial_number for atom in model.get_atoms()), default=0)
    #n_current_serial_number = max_serial_number + 1

    #get the chain
    chain = model[s_chain] #e.g. chain = model['A']
    l_residues = list(chain.get_residues())


    # Create and a new detached residue. we will set the specified position, but its not part of any chain yet
    new_residue = Residue.Residue((het_flag, n_residue, ' '), s_new_residue_name.ljust(3)[:3], '')


    for current_new_atom in d_new_residue_atoms:

        # define new atom in biopython
        atom = Atom.Atom(
            current_new_atom['structural_name'].upper().ljust(4),  # Atom name (e.g., 'C1', 'N1', 'O1') # Ensure atom name is 4 characters (padded or truncated)
            np.array(current_new_atom['coord'], dtype=float),  # Coordinates as a numpy array
            1.0,  # B-factor (optional, default 1.0)
            1.0,  # Occupancy (optional, default 1.0)
            ' ',  # Alternate location indicator
            str(current_new_atom['structural_name']).strip(),  # Atom full structural name. I dont know why this is necessary
            1001, #int(n_current_serial_number),  # serial number e.g. 1001 #i dont know why this is necessary, it seems like in the pdb file a renumbering takes place, and dont consider this
            str(current_new_atom['element_name']).strip().upper()# the element name (e.g., 'C', 'N', 'O')
        )
        #n_current_serial_number += 1
        new_residue.add(atom)
        #print(f"Atom: {atom.get_name()}, coord: {atom.get_coord()}, Residue: {s_new_residue_name}, Residueid: {n_residue}, Serial Number: {atom.serial_number}")

        
    # To insert the new resitue in the correct position, first detach all residues, saving all of them in a list
    residues_to_reinsert = []
    for residue in l_residues:
        residues_to_reinsert.append(residue)
        chain.detach_child(residue.id)

    #now recontruct the chain, putting the new residue in the correct position
    if s_replace_or_displace=='replace': #the new residue will replace the old one in that position

        for residue in residues_to_reinsert:
            if residue.id[1] == n_residue:
                chain.add(new_residue)#the new residue one is added
            else:
                chain.add(residue)#old residues remain the same

    elif s_replace_or_displace=='displace': #the new residue will be inserted in the position, displacing the folowing ones

        for residue in residues_to_reinsert:
            if residue.id[1] < n_residue:#current residue is before the new one
                chain.add(residue)#residues remain the same
            elif residue.id[1] == n_residue:#current residue is in the exact place of the new one
                chain.add(new_residue)#the new residue is added instead of the old one
            else:#current residue is after the new one
                previous_residue.id = (previous_residue.id[0], previous_residue.id[1] +1, previous_residue.id[2])#redefine the residue position, so the numbering will be correct in the pdb file
                chain.add(previous_residue)#residues remain the same, but is the one from the last iteration that will be inserted, so to displace all residues after the new one

                if residue.id[1] == len(residues_to_reinsert):#if it is the last residue, we inser also the final one and break the loop
                    residue.id = (residue.id[0], residue.id[1] +1, residue.id[2])#redefine the residue position, so the numbering will be correct in the pdb file
                    chain.add(residue)
                    break

            previous_residue = residue #save the current residue for the next iteration, if necessary



    # Write the modified structure to the output file
    io = PDBIO()
    io.set_structure(structure)
    io.save(s_output_pdb_file)
    print(f"Modified PDB saved to {s_output_pdb_file}")



def insert_atom_into_residue(s_input_pdb_file,s_output_pdb_file,s_chain,n_residue,s_atom_structural_name,s_atom_element_name,l_coord):
    """

    the user specify the chain and residue id. and put an atom there with defined name and coordinates

    
    cl.insert_atom_into_residue('pepticat6.pdb','pepticat6_with_atom.pdb','A',3,'CAE','AU',[1,1,1])

    """

    print("IN FUNCTION insert_atom_into_residue")

    
    


    # Parse the existing PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', s_input_pdb_file)

    # Get the first model or create a new one if needed
    model = structure[0] if len(structure) > 0 else print('CLEAN PIPE MESSAGE : ATTENTION there are several models in the pdb file, it is the first that will be edited')

    # Determine the maximum residue ID in the whole structure to avoid collisions
    max_serial_number = max((atom.serial_number for atom in model.get_atoms()), default=0)
    n_current_serial_number = max_serial_number + 1

    #get the chain
    chain = model[s_chain] #e.g. chain = model['A']
    l_residues = list(chain.get_residues())

    # Get the residue to be edited
    n_residue = n_residue-1#change residue number so it suits the list index
    residue = l_residues[n_residue] 


    # define new atom in biopython
    new_atom = Atom.Atom(
        s_atom_structural_name.upper().ljust(4),  # Atom name (e.g., 'C1', 'N1', 'O1') # Ensure atom name is 4 characters (padded or truncated)
        np.array(l_coord, dtype=float),  # Coordinates as a numpy array
        1.0,  # B-factor (optional, default 1.0)
        1.0,  # Occupancy (optional, default 1.0)
        ' ',  # Alternate location indicator
        str(s_atom_structural_name).strip(),  # Atom full structural name. I dont know why this is necessary
        1001, #int(n_current_serial_number),  # serial number e.g. 1001. I dont know why this is necessary. inside the new pdb, the numbering dont consider this
        str(s_atom_element_name).strip().upper() # the element name (e.g., 'C', 'N', 'O') it shoule be capital letters
    )


    # add atom to the residue
    residue.add(new_atom)


    # Write the modified structure to the output file
    io = PDBIO()
    io.set_structure(structure)
    io.save(s_output_pdb_file)
    print(f"Modified PDB saved to {s_output_pdb_file}")