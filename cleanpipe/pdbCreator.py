import PeptideBuilder
from Bio.PDB import PDBIO
import Geometry
from io import StringIO
from cleanpipe import topContent
from cleanpipe import filemanager
from cleanpipe import procedures
from cleanpipe import atomisticContruction
import subprocess
import os


def download_and_clean_pdb(s_molecule_name):
    """
    usage example:
    cl.download_and_clean_pdb("1aki")

    
    """

    #get pdb from portal
    subprocess.run(f"wget https://files.rcsb.org/download/{s_molecule_name}.pdb" , shell=True, check=True)

    #remove water
    subprocess.run(f"grep -v 'HOH' {s_molecule_name}.pdb > {s_molecule_name}_temp.pdb" , shell=True, check=True)
    subprocess.run(f"rm {s_molecule_name}.pdb" , shell=True, check=True)
    subprocess.run(f"mv {s_molecule_name}_temp.pdb {s_molecule_name}.pdb" , shell=True, check=True)


def create_peptide(s_outName, s_nTerminusCAP, s_aminoacids, s_cTerminusCAP, l_phi, l_psi_im1):

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
        atomisticContruction.add_acetyl_to_Nterminus(peptide)

    if s_cTerminusCAP == "amide":
        atomisticContruction.add_amide_to_Cterminus(peptide)


    #################################### create system. (ps this will add hydrogens) ###################################

    # Save temporary pdb file of the peptide
    io = PDBIO()
    io.set_structure(peptide)
    io.save(s_outName)




