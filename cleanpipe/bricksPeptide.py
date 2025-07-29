from cleanpipe import algelin
import Bio

import numpy as np
import PeptideBuilder
import Geometry

from Bio.PDB import Atom, Residue, Chain, PDBParser, PDBIO

from PeptideBuilder import Geometry
from PeptideBuilder import PeptideBuilder
import Bio.PDB



def add_ACE_to_Nterminus(peptide):
    
    chain = peptide[0]['A']
    l_residues = list(chain.get_residues())

    # Get the first residue in the chain and its coordinates
    first_residue = l_residues[0]
    coordsN = first_residue['N'].get_coord()
    coordsCA = first_residue['CA'].get_coord()
    coordsC = first_residue['C'].get_coord()

    # Create the first C of ACE
    coords = algelin.find_new_atom_coord(coordsN, coordsCA, coordsC, 1.5, -60, 0)
    acetyl_c1 = Bio.PDB.Atom.Atom("C", coords, 0.0, 1.0, ' ', 'C', 1001, 'C')

    # Create the second C, and the O of ACE
    coords = algelin.find_new_atom_coord(acetyl_c1.coord, coordsN, coordsCA, 1.5, -45, -90)
    acetyl_c2 = Bio.PDB.Atom.Atom("CH3", coords, 0.0, 1.0, ' ', 'CH3', 1002, 'C')
    coords = algelin.find_new_atom_coord(acetyl_c1.coord, coordsN, coordsCA, 1.5, 75, 0)
    acetyl_o = Bio.PDB.Atom.Atom("O", coords, 0.0, 1.0, ' ', 'O', 1003, 'O')

    #IMPORTANT: THE NAMES OF THE ATOMS IN ACE HAVE TO MATCH THE NAMES OF THE FORCEFIELD YOU WILL CHOSE IN THE FUTURE
    #THE NAMES "CH3", "C" AND "CA" COMM FROM "CHARMM36". THEY CAN BE FOUND IN THE FILE "aminoacids.hdb" 

    # Create the new ACE residue
    ace_residue = Bio.PDB.Residue.Residue((' ', 1, ' '), 'ACE', '    ')
    ace_residue.add(acetyl_c1)
    ace_residue.add(acetyl_c2)
    ace_residue.add(acetyl_o)

    # Detach all residues to avoid index conflict
    for residue in l_residues:
        chain.detach_child(residue.id)

    # Insert the ACE residue
    chain.add(ace_residue)

    # Re-insert each residue with updated numbers
    for i, residue in enumerate(l_residues, start=2):
        residue.id = (residue.id[0], i, residue.id[2])
        chain.add(residue)

def add_NME_to_Cterminus(peptide):
    
    chain = peptide[0]['A']
    l_residues = list(chain.get_residues())
    

    # Get the first residue in the chain and its coordinates
    last_residue = l_residues[len(l_residues)-1]
    coordsC = last_residue['C'].get_coord()
    coordsCA = last_residue['CA'].get_coord()
    coords0 = last_residue['O'].get_coord()

    # Create the N of NME
    coords = algelin.find_new_atom_coord(coordsC, coordsCA, coords0, 1.5, -45, 90)
    acetyl_n = Bio.PDB.Atom.Atom("N", coords, 0.0, 1.0, ' ', 'N', 1001, 'N')

    # Create the C of NME
    coords = algelin.find_new_atom_coord(acetyl_n.coord, coordsC, coords0, 1.5, +45, 0)
    acetyl_c = Bio.PDB.Atom.Atom("CH3", coords, 0.0, 1.0, ' ', 'CH3', 1002, 'C')


    #IMPORTANT: THE NAMES OF THE ATOMS IN ACE HAVE TO MATCH THE NAMES OF THE FORCEFIELD YOU WILL CHOSE IN THE FUTURE
    #THE NAMES "CH3", "C" AND "O" COMM FROM "CHARMM36". THEY CAN BE FOUND IN THE FILE "aminoacids.hdb" 

    # Create the new NME residue
    nme_residue = Bio.PDB.Residue.Residue((' ', len(l_residues)+1, ' '), 'NME', '    ')
    nme_residue.add(acetyl_n)
    nme_residue.add(acetyl_c)


    # Insert the NME residue
    chain.add(nme_residue)









def add_truss(peptide, p1, p2):
    """

    xxx there is a duplicate function in bricksPDB

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





    #add all that truss to the fisrt residue of the peptide
    chain = peptide[0]['A']
    l_residues = list(chain.get_residues())
    residue1 = l_residues[1]

    idvert = 1
    for vertice in vertices:
        residue1.add(Bio.PDB.Atom.Atom("X"+str(idvert), vertice, 0.0, 1.0, ' ', "X"+str(idvert), 100+idvert, ''))
        idvert += 1


    #in addition to the inplace modification of the peptide object, I will return the vertices and the distances
    return vertices, distances
