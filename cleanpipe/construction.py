import PeptideBuilder
from Bio.PDB import PDBIO
import Geometry
from io import StringIO
from cleanpipe import construction_aux
from cleanpipe import topContent
from cleanpipe import filemanager
import subprocess
import os







def create_peptide_in_solution(s_outName, s_nTerminusCAP, s_aminoacids, s_cTerminusCAP, l_phi, l_psi_im1, s_solvent, s_forceField, s_boxSize):
    """
    example:
    create_peptide_in_solution("alaHW","acyl","AAAAAA","amide",[-57.8,-57.8,-57.8,-57.8,-57.8,-57.8],[-47.0,-47.0,-47.0,-47.0,-47.0,-47.0],"tip3p","charmm36-jul2022", "5.1 5.1 5.1")
    
    this are the N terminus CAP 
    this are the C terminus CAP 

    s_outFileName   : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    s_nTerminusCAP  : options: "acyl" , none
    s_aminoacids    : string of aminoacid one letter code, for example "AAAAGGAALL"
    s_cTerminusCAP  : options: "amide", none
    l_phi           : vector with angles
    l_psi_im1       : vector with angles
    s_solvent       : choose a water model, such as "tip3p", or put a folder containg a system with a solvent box, sucha as "octn_filledbox", making sure there are a octn_filledbox.gro and octn_filledbox.top inside that folder
    s_forceField    : one of the gromacs recognized force fields
    s_boxSize       : string with x y z sizes, for example "3 3 3"
     
    """

    subprocess.run(f"mkdir {s_outName}", shell=True)
    s_outPathAndName = f"{s_outName}/{s_outName}"
    



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


    #################################### create system. this will add hydrogens ###################################

    # Save temporary pdb file
    io = PDBIO()
    io.set_structure(peptide)
    io.save(f"{s_outPathAndName}.pdb")

    #pdb2gmx
    subprocess.run(f"printf '8\n7\n' | gmx pdb2gmx -f {s_outPathAndName}.pdb -o {s_outPathAndName}.gro -p {s_outPathAndName}.top -i {s_outPathAndName}.posres.itp -missing -ter -ignh -water tip3p -ff {s_forceField}", shell=True)
    
    #pdb is not necessary anymore. as pdb2gmx provided the gro and top files to describe the system
    subprocess.run(f"rm {s_outPathAndName}.pdb" , shell=True)# I chose to overwrite the old gro
    
    #define box size. s_boxSize contains the user definition (ex: "3 3 3")
    subprocess.run(f"gmx editconf -f {s_outPathAndName}.gro -o {s_outPathAndName}.gro -c -box {s_boxSize} -bt cubic", shell=True)
    subprocess.run(f"rm {s_outName}/\\#{s_outName}.gro.1\\#" , shell=True)# I chose to overwrite the old gro


    ################################## add solvent to the system. I have 2 options here: tip3p or filled box ##############################################


    if s_solvent == "tip3p":
        #this mean the user has chosen the tip3p water model, already part of the forcefield

        subprocess.run(f"gmx solvate -cp {s_outPathAndName}.gro -p {s_outPathAndName}.top -o {s_outPathAndName}.gro", shell=True)
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.gro.1\\#" , shell=True)#I chose to overwrite the old gro
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.top.1\\#" , shell=True)#I chose to overwrite the old top

        #set the the name of the system in the top file
        topContent.setSystemName(f"{s_outPathAndName}.top", f"peptide in water (tip3p)" )


    elif filemanager.check_folder_solvent_box(s_solvent) == True:
        #this mean the user has chosen a folder (ex: path/to/folder/octn) containing a system that is a box filled with solvent. for example octn.gro and octn.itp

        #get the base name (ex: octn)
        s_sol_base_name = os.path.basename(os.path.normpath(s_solvent))

        #edit gro to insert the solvent
        subprocess.run(f"gmx solvate -cp {s_outPathAndName}.gro -cs {s_solvent}//{s_sol_base_name}.gro -p {s_outPathAndName}.top -o {s_outPathAndName}.gro", shell=True)
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.gro.1\\#" , shell=True)#I chose to overwrite the old gro
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.top.1\\#" , shell=True)#I chose to overwrite the old top

        #edit top to insert a line including a reference of the solvent itp before the [ system ] directive
        subprocess.run(rf'''awk -v line='#include "{s_sol_base_name}.itp"' '/\[ system \]/{{print line"\n"; i=2}}i&&!--i{{next}}1' {s_outPathAndName}.top > temp.top && mv temp.top {s_outPathAndName}.top''', shell=True, check=True)

        #make that itp file realy available in the current system
        subprocess.run(f"cp {s_solvent}/{s_sol_base_name}.itp {s_outName}/" , shell=True)

        #set the the name of the system in the top file
        topContent.setSystemName(f"{s_outPathAndName}.top", f"peptide in custom solvent" )

    else:
        print("solvation failed")
