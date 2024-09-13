import PeptideBuilder
from Bio.PDB import PDBIO
import Geometry
from io import StringIO
from cleanpipe import construction_aux
from cleanpipe import topContent
from cleanpipe import filemanager
from cleanpipe import procedures
import subprocess
import os



def create_peptide_in_solution(s_outName, s_nTerminusCAP, s_aminoacids, s_cTerminusCAP, l_phi, l_psi_im1, s_solvent, s_forceField, s_boxSize):
    """
    usage example:
    cl.create_peptide_in_solution("ala_in_solvent","acyl","AAAAAA","amide",[-57.8,-57.8,-57.8,-57.8,-57.8,-57.8],[-47.0,-47.0,-47.0,-47.0,-47.0,-47.0],"tip3p","charmm36-jul2022", "5.1 5.1 5.1")
    
    will create the peptide and the entire system out of nowere (no input files required)
    a folder containing the system will be created. this will be done using just the function arguments 

    s_outFileName   : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    s_nTerminusCAP  : N terminus CAP. there are only the folowing options: "acyl" , none
    s_aminoacids    : string of aminoacid one letter code, for example "AAAAGGAALL"
    s_cTerminusCAP  : this are the C terminus CAP. there are only the folowing options:: "amide", none
    l_phi           : vector with angles
    l_psi_im1       : vector with angles
    s_solvent       : choose a water model, for example as "tip3p", or a folder, for example "octn_filledbox". The folder have to contain a system with a solvent box, in other words, it has to contain a octn_filledbox.gro and a octn.itp
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"
     
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
        construction_aux.add_acetyl_to_Nterminus(peptide)

    if s_cTerminusCAP == "amide":
        construction_aux.add_amide_to_Cterminus(peptide)


    #################################### create system. (ps this will add hydrogens) ###################################

    # Save temporary pdb file of the peptide
    io = PDBIO()
    io.set_structure(peptide)
    io.save("temp.pdb")

    # Call this important function, that will construct the entire system aroud the peptide
    procedures.pdb2molecule_in_solvent("temp.pdb", s_outName, s_solvent, s_forceField, s_boxSize)

    # temporary pdb of the peptide is not necessary anymore
    subprocess.run(f"rm temp.pdb" , shell=True)


    ################################### set the the name of the system in the top file #######################################
    s_outPathAndName = f"{s_outName}/{s_outName}"
    topContent.setSystemName(f"{s_outPathAndName}.top", f"{s_outName} (custom peptide, insterted in solution made using {s_solvent})" )
