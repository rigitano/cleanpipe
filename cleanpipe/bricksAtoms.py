from cleanpipe import algelin
import Bio



def add_acetyl_to_Nterminus(peptide):
    
    chain = peptide[0]['A']
    l_residues = list(chain.get_residues())

    # Get the first residue in the chain and its coordinates
    first_residue = l_residues[0]
    coordsN = first_residue['N'].get_coord()
    coordsCA = first_residue['CA'].get_coord()
    coordsCB = first_residue['CB'].get_coord()

    # Create the first C of ACE
    coords = algelin.find_new_atom_coord(coordsN, coordsCA, coordsCB, 1.5, -60, 0)
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

def add_amide_to_Cterminus(peptide):
    
    chain = peptide[0]['A']
    l_residues = list(chain.get_residues())
    

    # Get the first residue in the chain and its coordinates
    last_residue = l_residues[len(l_residues)-1]
    coordsC = last_residue['C'].get_coord()
    coordsCA = last_residue['CA'].get_coord()
    coords0 = last_residue['O'].get_coord()

    # Create the first C of ACE
    coords = algelin.find_new_atom_coord(coordsC, coordsCA, coords0, 1.5, -45, 90)
    acetyl_n = Bio.PDB.Atom.Atom("N", coords, 0.0, 1.0, ' ', 'N', 1001, 'N')

    # Create the second C, and the O of ACE
    coords = algelin.find_new_atom_coord(acetyl_n.coord, coordsC, coords0, 1.5, +45, 0)
    acetyl_c = Bio.PDB.Atom.Atom("CH3", coords, 0.0, 1.0, ' ', 'CH3', 1002, 'C')


    #IMPORTANT: THE NAMES OF THE ATOMS IN ACE HAVE TO MATCH THE NAMES OF THE FORCEFIELD YOU WILL CHOSE IN THE FUTURE
    #THE NAMES "CH3", "C" AND "O" COMM FROM "CHARMM36". THEY CAN BE FOUND IN THE FILE "aminoacids.hdb" 

    # Create the new ACE residue
    nme_residue = Bio.PDB.Residue.Residue((' ', len(l_residues)+1, ' '), 'NME', '    ')
    nme_residue.add(acetyl_n)
    nme_residue.add(acetyl_c)


    # Insert the NME residue
    chain.add(nme_residue)