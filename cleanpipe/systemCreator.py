from cleanpipe import pdbCreator
from cleanpipe import filemanager
from cleanpipe import topContent
from cleanpipe import betterGromacs

import subprocess
import re
import os



def pdb2filled_box(s_pdbfile, s_forceField):
    """

    usage example:
    cl.pdb2filled_box("octn.pdb","charmm36-jul2022")

    create a 5x5x5 box system filled with a lot of copies of the molecule

    as gromacs dont have such tool, its necessary to be creative, and use "pdb2gmx" to create a system with 1 molecule, then
    use "insert-molecules" to fill the box with copyes of the molecule, then modify the top to reflect the new total
    the top file is also edited to change name of the system. and also the name of the molecule
    here we also make sure the outputs are renamed to be blabla_filledbox.gro, blabla_filledbox.top and blabla_filledbox.posres.top

    """

    #check if the filename inside s_pdbfile is valid
    filemanager.check_file(s_pdbfile,['.pdb']) 
    #obtein just the file name. ex: blabla/blabla/filename.bla
    s_filename = filemanager.get_filename_without_extension(s_pdbfile) 


    subprocess.run(f"mkdir {s_filename}_filled_box", shell=True, check=True)
    s_outPathAndName = f"{s_filename}_filled_box/{s_filename}"


    #create a system with 1 molecule.
    subprocess.run(f"gmx pdb2gmx -f {s_filename}.pdb -o {s_outPathAndName}.gro -p {s_outPathAndName}.top -i posres.itp -water none -ff {s_forceField}" , shell=True, check=True)
    
    #pdb2gmx generates a useless posres.itp with useless posres for 1 molecule. so I delete the posres.itp and the inclusion in the top
    subprocess.run(f"rm posres.itp" , shell=True, check=True) 
    topContent.remove_posres_inclusion(f"{s_outPathAndName}.top")

    #manipulate the GRO file to create a 5x5x5 box and fill it with copyes of the molecule
    result = subprocess.run(f"gmx insert-molecules -ci {s_outPathAndName}.gro -nmol 1000 -rot -box 5 5 5 -o {s_outPathAndName}_filled_box.gro" , shell=True, check=True, capture_output=True,text=True)# 
    print(result.stdout+result.stderr)
    subprocess.run(f"rm {s_outPathAndName}.gro" , shell=True, check=True)# now that we have the filled box gro, the 1 molecule gro can be deleted
    print(f"gro file written: \n                      {s_outPathAndName}_filled_box.gro")

    #get the number of added molecules. 
    match = re.search(r'Added\s+(\d+)\s+molecules', result.stdout+result.stderr)
    added_molecules = int(match.group(1))

    #rename the top and posres.itp files.
    subprocess.run(f"mv {s_outPathAndName}.top {s_outPathAndName}_filled_box.top" , shell=True, check=True)

    #change the ugly molecule name currently inside the TOP file.
    uglyMolName = topContent.getMoleculeName(f"{s_outPathAndName}_filled_box.top")
    molName = s_filename
    topContent.replaceMoleculeName(f"{s_outPathAndName}_filled_box.top", uglyMolName, molName)

    #update the TOP file with the new total the molecule
    topContent.update_molecule_quantity(f"{s_outPathAndName}_filled_box.top", molName, added_molecules)

    #split the TOP file, into a ITP that describes the molecule and a simple TOP that contains only name of the system and the totals.
    topContent.decompose_TOP_file_into_SOCKETTOP_and_ITPs(f"{s_outPathAndName}_filled_box.top")
    

    #give a name for the system
    topContent.setSystemName(f"{s_outPathAndName}_filled_box.top", f"box filled with {s_filename}" )


def pdb2molecule_in_solvent(s_pdbfile, s_outSytemName, s_solvent, s_forceField, s_boxSize):
    """
    usage example:
    cl.pdb2molecule_in_solvent("1LZ1.pdb", "1LZ1_in_water", "tip3p", "charmm36-jul2022", "3 3 3")
    cl.pdb2molecule_in_solvent("1LZ1.pdb", "1LZ1_in_octane", "octn_filled_box", "charmm36-jul2022", "3 3 3")
    

    s_pdbfile       : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
    s_outSytemName  : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    s_solvent       : choose a water model, for example as "tip3p", or a folder, for example "octn_filledbox". The folder have to contain a system with a solvent box, in other words, it has to contain a octn_filledbox.gro and a octn.itp
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"
    """

    # check if the filename inside s_pdbfile is valid
    filemanager.check_file(s_pdbfile,['.pdb']) 

    # create gro and top from pdb. then add the box size to the gro
    betterGromacs.better_pdb2gmx(s_pdbfile,s_outSytemName,s_forceField,s_boxSize)

    # add solvent to the system. I have 2 options here: tip3p or filled box
    betterGromacs.better_solvate(s_outSytemName,s_solvent)

    # set the the name of the system in the top file 
    topContent.setSystemName(f"{s_outSytemName}/{s_outSytemName}.top", f"{s_outSytemName} (molecule from {s_pdbfile}, inserted in solution made using {s_solvent})" )




def void2peptide_in_solvent(s_peptideName, s_systemName, s_nTerminusCAP, s_aminoacids, s_cTerminusCAP, l_phi, l_psi_im1, s_solvent, s_forceField, s_boxSize):
    """
    usage example:
    cl.void2peptide_in_solvent("poliA","poliA_in_water","acyl","AAAAAA","amide",[-57.8,-57.8,-57.8,-57.8,-57.8,-57.8],[-47.0,-47.0,-47.0,-47.0,-47.0,-47.0],"tip3p","charmm36-jul2022", "5.1 5.1 5.1")
    cl.void2peptide_in_solvent("poliA","poliA_in_octane","acyl","AAAAAA","amide",[-57.8,-57.8,-57.8,-57.8,-57.8,-57.8],[-47.0,-47.0,-47.0,-47.0,-47.0,-47.0],"octn_filled_box","charmm36-jul2022", "5.1 5.1 5.1")


    will create the peptide and the entire system out of nowere (no input files required)
    a folder containing the system will be created. this will be done using just the function arguments 

    s_peptideName   : name of the peptide that will be created, for example alaH
    s_systemName    : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    s_nTerminusCAP  : N terminus CAP. there are only the folowing options: "acyl" , none
    s_aminoacids    : string of aminoacid one letter code, for example "AAAAGGAALL"
    s_cTerminusCAP  : this are the C terminus CAP. there are only the folowing options:: "amide", none
    l_phi           : vector with angles
    l_psi_im1       : vector with angles
    s_solvent       : choose a water model, for example as "tip3p", or a folder, for example "octn_filledbox". The folder have to contain a system with a solvent box, in other words, it has to contain a octn_filledbox.gro and a octn.itp
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"
     
    """

    # create  pdb file containing only a peptide. in this case the pdb will be temporary, and deleted after the conversion to top and gro
    pdbCreator.create_peptide(f"{s_peptideName}.pdb", s_nTerminusCAP, s_aminoacids, s_cTerminusCAP, l_phi, l_psi_im1)
   
    # create gro and top from pdb. then add the box size to the gro. this is the my improved gromacs version that create the system in a new folder
    betterGromacs.better_pdb2gmx(f"{s_peptideName}.pdb",s_systemName,s_forceField,s_boxSize, b_addterminal = False)# this last parameter was set to False so not to add termini, as they are already present in the peptide in this case)

    # temporary pdb of the peptide is not necessary anymore
    subprocess.run(f"rm {s_peptideName}.pdb" , shell=True, check=True)

    # add solvent to the system. this is my improved gromacs, the imput is a folderthere are 2 possible options here: tip3p or filled box
    betterGromacs.better_solvate(s_systemName,s_solvent)

    # set the the name of the system in the top file 
    topContent.setSystemName(f"{s_systemName}/{s_systemName}.top", f"{s_systemName} (custom peptide, insterted in solution made using {s_solvent})" )