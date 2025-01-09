import PeptideBuilder
from Bio.PDB import PDBIO
import Geometry
from cleanpipe import bricksFileSystem
from cleanpipe import bricksAtoms
import subprocess

# xxx change name to bricksPdb


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


def create_peptide(s_outName, s_aminoacids, l_phi, l_psi_im1, s_nTerminusCAP, s_cTerminusCAP, ):

    """
    
    example:
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
    if s_nTerminusCAP == "acyl":
        bricksAtoms.add_acetyl_to_Nterminus(peptide)

    if s_cTerminusCAP == "amide":
        bricksAtoms.add_amide_to_Cterminus(peptide)


    #################################### create system. (ps this will add hydrogens) ###################################

    # Save temporary pdb file of the peptide
    io = PDBIO()
    io.set_structure(peptide)
    io.save(s_outName)




from Bio.PDB import PDBParser, PDBIO, Atom, Residue, Chain, Model, Structure
import numpy as np

def insert_new_molecule_into_pdb(input_pdb_file, output_pdb_file, new_molecule_atoms, new_chain_id='Z'):
    """
    Insert a new molecule into a pre-existing PDB file as a new chain.

    Parameters:
        input_pdb_file (str): Path to the input PDB file.
        output_pdb_file (str): Path to save the modified PDB file.
        new_molecule_atoms (list of dict): Each dict contains 'name', 'position', and 'residue_name'.
        new_chain_id (str): The ID of the new chain to create.



    # Example 
    new_molecule_atoms = [
        {'name': 'C1', 'position': [15.0, 12.3, 10.7], 'residue_name': 'X01'},
        {'name': 'N1', 'position': [15.5, 13.0, 11.2], 'residue_name': 'X01'},
        {'name': 'O1', 'position': [16.0, 14.0, 12.1], 'residue_name': 'X02'},
        {'name': 'H1', 'position': [14.5, 11.5, 9.9], 'residue_name': 'H2O'}
    ]

    cl.insert_new_molecule_into_pdb("example.pdb", "modified_example.pdb", new_molecule_atoms) 
    """

    print("IN FUNCTION insert_new_molecule_into_pdb")
    print(f"Inserting new molecule into PDB file {input_pdb_file} as chain {new_chain_id}")
    print(new_molecule_atoms)

    # Parse the existing PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', input_pdb_file)

    # Get the first model or create a new one if needed
    model = structure[0] if len(structure) > 0 else Model.Model(0)

    # Add a new chain for the molecule
    new_chain = Chain.Chain(new_chain_id)
    model.add(new_chain)

    # Determine the maximum residue ID in the whole structure to avoid collisions
    max_residue_id = max((residue.id[1] for residue in model.get_residues()), default=0)
    max_serial_number = max((atom.serial_number for atom in model.get_atoms()), default=0)

    # Insert residues and atoms into the new chain
    current_residue_id = max_residue_id + 1
    current_serial_number = max_serial_number + 1


    for atom_info in new_molecule_atoms:
        residue_name = atom_info['residue_name']

        # Automatically assign residue IDs sequentially
        residue_id = (' ', current_residue_id, ' ')  # (het_flag, resseq, icode)

        # Create and add a new residue
        residue = Residue.Residue(residue_id, residue_name, ' ')
        new_chain.add(residue)

        # Create new atom
        atom = Atom.Atom(
            atom_info['name'],  # Atom name (e.g., 'C1', 'N1', 'O1')
            tuple(atom_info['position']),  # Coordinates as a tuple (x, y, z)
            1.0,  # B-factor (optional, default 1.0)
            1.0,  # Occupancy (optional, default 1.0)
            '',  # Alternate location indicator
            atom_info['name'],  # Full atom name
            element=atom_info['name'][0],  # Guess element from the first character of the atom name
            serial_number=current_serial_number  # Automatically assigned serial number
        )

        # Add atom to the residue
        residue.add(atom)

        # Increment residue ID for the next residue
        current_residue_id += 1
        current_serial_number += 1

    # Write the modified structure to the output file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file)
    print(f"Modified PDB saved to {output_pdb_file}")



def add_truss2(s_pdb_file, p1, p2):
    """

    vertices, edges = add_truss(peptide,[0, 0, 0], [1, 1, 1])
    print(vertices)
    print(edges)
    """

    #define the number of equidistant points between the atips, and the size of the square
    n = 6
    square_size = 6

    # Calculate the direction vector and the step size
    direction = np.array(p2, dtype=float) - np.array(p1, dtype=float)
    step = direction / (n + 1)
    
    # Normalize the direction vector
    direction = direction / np.linalg.norm(direction)
    
    # Generate the equidistant points
    points = [np.array(p1, dtype=float) + i * step for i in range(n + 2)]
    
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
            vertex = point + (dx * square_size / 2) * perp_vector1 + (dy * square_size / 2) * perp_vector2
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

    #inset the coordinates in a pdb file
    new_molecule_atoms = []
    for i in range(len(vertices)):
        new_molecule_atoms.append({'name': f"X{i}", 'position': vertices[i], 'residue_name': 'TRS'})#example: {'name': 'X1', 'position': [15.0, 12.3, 10.7], 'residue_name': 'TRS'}
        print(new_molecule_atoms)

    #at last, we insert the new molecule in the pdb
    insert_new_molecule_into_pdb(s_pdb_file, "hahaha.pdb", new_molecule_atoms, new_chain_id='T')

