from cleanpipe import bricksPDB
from cleanpipe import bricksFileSystem
from cleanpipe import bricksTOP
from cleanpipe import bricksGROMACS

import subprocess
import re
import sys
import os
import functools
import shutil

from pathlib import Path


def ensure_original_directory(func):
    """
    this is a function to use as a decorator in all the functions that change the directory where the python process is run.
    this happens, for example in some functions that use gromacs but I want the output to be saved in a new folder, created just beside the input file. 
    this decorator guarantees we go out of that folder even if the function crashes
    
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Store the original directory
        original_directory = os.getcwd()
        try:
            # Execute the function
            return func(*args, **kwargs)
        finally:
            # Return to the original directory, even if an error occurred
            os.chdir(original_directory)
    return wrapper

def pdb2box_full_of_that(s_pdbfile, s_forceField, s_box_size, n_mol_max):
    """

    usage example:
    cl.pdb2box_full_of_that("octn.pdb","charmm36-jul2022", "5 5 5", 1000)

    create a 5x5x5 box system filled with a lot of copies of the molecule

    as gromacs dont have such tool, its necessary to be creative, and use "pdb2gmx" to create a system with 1 molecule, then
    use "insert-molecules" to fill the box with copyes of the molecule, then modify the top to reflect the new total
    the top file is also edited to change name of the system. and also the name of the molecule
    here we also make sure the outputs are renamed to be blabla_filledbox.gro, blabla_filledbox.top and blabla_filledbox.posres.top

    """

    #check if the filename inside s_pdbfile is valid
    bricksFileSystem.check_extention(s_pdbfile,['.pdb']) 
    #obtain just the file name. ex: blabla/blabla/filename.bla
    s_filename = bricksFileSystem.get_filename_without_extension(s_pdbfile) 


    bricksFileSystem.run_and_capture(f"mkdir box_full_of_{s_filename}")
    s_outPathAndName = f"box_full_of_{s_filename}/box_full_of_{s_filename}"


    #create a system with 1 molecule.
    bricksFileSystem.run_and_capture(f"gmx pdb2gmx -f {s_filename}.pdb -o {s_outPathAndName}_just1mol.gro -p {s_outPathAndName}.top -i posres.itp -water none -ff {s_forceField}")
    
    #pdb2gmx generates a useless posres.itp with useless posres for 1 molecule. so I delete the posres.itp and the inclusion in the top
    bricksFileSystem.delete("posres.itp")
    bricksTOP.remove_posres_inclusion(f"{s_outPathAndName}.top")

    #manipulate the GRO file to create a and fill it with copyes of the molecule
    captured_output = bricksFileSystem.run_and_capture(f"gmx insert-molecules -ci {s_outPathAndName}_just1mol.gro -nmol {str(n_mol_max)} -rot xyz -box {s_box_size} -o {s_outPathAndName}.gro")
    print(f"\nCLEANPIPE MESSAGE\ngro file written: \n                     {s_outPathAndName}.gro")

    #now we have the final gro with a lot of molecules. its time to delete the initial one
    bricksFileSystem.delete(f"{s_outPathAndName}_just1mol.gro")

    #get the number of added molecules. 
    match = re.search(r'Added\s+(\d+)\s+molecules', captured_output)
    added_molecules = int(match.group(1))

    #change the ugly molecule name currently inside the TOP file.
    uglyMolName = bricksTOP.getMoleculeName(f"{s_outPathAndName}.top")
    molName = s_filename
    bricksTOP.replaceMoleculeName(f"{s_outPathAndName}.top", uglyMolName, molName)

    #update the TOP file with the new total the molecule
    bricksTOP.update_molecule_quantity(f"{s_outPathAndName}.top", molName, added_molecules)

    #split the TOP file, into a ITP that describes the molecule and a simple TOP that contains only name of the system and the totals.
    bricksTOP.decompose_TOP_file_into_TOP_and_ITPs(f"{s_outPathAndName}.top")
    

    #give a name for the system
    bricksTOP.setSystemName(f"{s_outPathAndName}.top", f"box filled with {s_filename}" )

@ensure_original_directory
def pdb2molecule_in_solvent(s_pdbfile, s_outSytemName, solvent, s_forceField, s_boxSize, s_maxsol=0, s_aditional_arguments=''):
    """
    s_pdbfile       : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
    s_outSytemName  : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    solvent         : choose a water model, for example as "tip3p", or a list containg a gro of a box of solvents and its itp, for example ["../solvents/box_full_of_octn.gro","octn.itp], bot the file name should be set in reference to the system top
    s_forceField    : place where the ff is, relative to the current location, so to be copyed to the system folder. examples:
                                                                                                                            '../alanine12/v15-truss/charmm36-jul2022.ff'
                                                                                                                            '~/ff/martini3001'
    s_boxSize       : string with x y z sizes, for example "3 3 3"
    s_maxsol        : the maximum number of solvent molecules that will be added. this is optional here, as the 0 value mean the parameter wont be considered by gromacs

    s_aditional_arguments : add extra arguments in martinize or pdb2gmx. for example, in martinize: "-water-bias -water-bias-eps E:-0.5 H:-1.0 -ss HHHHHHHHHHHH"

    examples:

    cl.pdb2molecule_in_solvent("g12H.pdb", "g12HW", "tip3p", 'charmm36-jul2022', "6 6 6", "6943")
    cl.pdb2molecule_in_solvent("g12H.pdb", "g12HO", ["../solvents/box_full_of_octn.gro","../solvents/octn.itp], 'charmm36-jul2022', "6.1 6.1 6.1", "712")

    cl.pdb2molecule_in_solvent("g12H.pdb", "g12HW_cg", ["~/solvents/martini/Water-pure/water.gro","~/ff/martini3001/martini_v3.0.0_solvents_v1.itp"], '~/ff/martini3001', "6 6 6", "1736", "-ss HHHHHHHHHHHH")
    cl.pdb2molecule_in_solvent("g12H.pdb", "g12HO_cg", ["~/solvents/martini/Octane/OCT_PRO1.gro", "~/ff/martini3001/martini_v3.0.0_solvents_v1.itp"], '~/ff/martini3001', "6.1 6.1 6.1", "712", "-ss HHHHHHHHHHHH")

    """
    print("CLEANPIPE called function pdb2molecule_in_solvent")

    # check if the filename inside s_pdbfile is valid
    bricksFileSystem.check_extention(s_pdbfile,['.pdb']) 


    ########################### create gro, itps and top in a folder with the name of the system ###########################
    bricksGROMACS.pdb2system(s_pdbfile, s_outSytemName, s_forceField, s_boxSize, s_aditional_arguments)
    ########################################################################################################################


    ######################### add solvent to the system. I have 2 options here: tip3p or filled box ########################
    bricksGROMACS.solvate_and_neutralize(s_outSytemName, solvent, s_maxsol)
    ########################################################################################################################


    # set the the name of the system in the top file 
    if isinstance(solvent, str): #the user inserted a string, that should mean a water model (ex "tip3p")
        s_solvent_text = solvent + ".gro"
    elif isinstance(solvent, list): #the user inserted a list, that should mean a solvent gro and itps (ex ["../solvents/box_full_of_octn.gro","octn.itp])
        s_solvent_text = bricksFileSystem.get_filename_with_extension(solvent[0])
    else:
        print("CLEANPIPE error, solvent content is unexpected")
    bricksTOP.setSystemName(f"{s_outSytemName}/{s_outSytemName}.top", f"{s_outSytemName} ; molecule from \"{s_pdbfile}\", inserted in solvent from \"{s_solvent_text}\"" )


    #copy usefull scripts to the system folder
    files_to_copy = [
        "runFEPoff.sh",
        "runREALISTIC.sh",
    ]
    module_path = Path(__file__).resolve().parent # Where the cl module lives
    source_dir = module_path / "bash" # Folder containing the source files
    dest_dir = Path(s_outSytemName) # Destination folder
    for filename in files_to_copy:
        src = source_dir / filename
        dst = dest_dir / filename
        shutil.copy(src, dst)




@ensure_original_directory
def pdb2molecule_in_water_and_octane(s_pdbfile, s_folderName, s_forceField, s_boxSize, s_maxsolW, s_maxsolO, s_aditional_arguments=''):
    """
    this function is perfect for canclulations of free energies of transfer. there will be two systems in the same folder. 
    but for DG transfer, this actualy makes things easyer. you can run two runREALISTIC scrits, and then two runFEPoff scripts. and all the data you need will be there

    s_pdbfile       : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
    s_folderName    : string with the name of the folder, for example "ala6Hdihr200-transfer". this folder will be created, and inside it, the two systems will be created inside it
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"
    s_maxsolW       : max number of water molecules. ex 
                                                           1residues  and box 3 3 3 : 873 (xxx for martini)
                                                           6residues  and box 5 5 5 : 3616?
                                                           11residues and box 6 6 6 : 6943
                                                           12residues and box 6 6 6 : 6943 (1736 for martini)
                                                           13residues and box 6 6 6 : xxx (xxx for martini)

    s_maxsolO       : max number of octane molecules. ex:
                                                           1residues  and box 3 3 3 : 65
                                                           6residues  and box 5 5 5 : 395?
                                                           11residues and box 6 6 6 : 712
                                                           12residues and box 6 6 6 : 712
                                                           13residues and box 6 6 6 : xxx (xxx for martini)

    s_aditional_arguments : extra stuff you might want to add to martinize2 or pdb2gmx

    example:
    cl.pdb2molecule_in_water_and_octane("normal_peptide.pdb", "ala6Hdih200-transfer", "charmm36-jul2022", "6 6 6", "6943", "712" )
    cl.pdb2molecule_in_water_and_octane("g12H.pdb",           "g12Hbias1-transfer",   '~/ff/martini3001', "6 6 6", "1736", "712", s_aditional_arguments='-water-bias -water-bias-eps E:-0.5 H:-1.0 -ss HHHHHHHHHHHH')


    """
    print("CLEANPIPE called function pdb2molecule_in_water_and_octane")
    module_path = Path(__file__).resolve().parent # Where the cl module lives
    ff_dir       = module_path / "USEFUL_FORCEFIELDS" 
    solvents_dir = module_path / "USEFUL_SOLVENTS" 

    if "martini3001" in s_forceField.lower():
        solvent1  = [ solvents_dir / "martini3001" / "Water-pure" / "water.gro",     ff_dir / "martini3001" / "martini_v3.0.0_solvents_v1.itp" ]    # I think water is already in the ff itp, so no need for this solvents file       
        solvent2  = [ solvents_dir / "martini3001" / "Octane"     / "OCT_PRO1.gro",  ff_dir / "martini3001" / "martini_v3.0.0_solvents_v1.itp" ]
        s_boxSize1 = s_boxSize #"6 6 6"
        s_boxSize2 = s_boxSize #"6.1 6.1 6.1"
        s_maxsol1 = s_maxsolW
        s_maxsol2 = s_maxsolO
    elif "martini22" in s_forceField.lower():
        solvent1  = [ solvents_dir / "martini22" / "Water-pure" / "water.gro",     ff_dir / "martini22" / "martini_v2.0_solvents.itp" ] # I think water is already in the ff itp, so no need for this solvents file  
        solvent2  = [ solvents_dir / "martini22" / "Octane"     / "OCT_PRO1.gro",  ff_dir / "martini22" / "martini_v2.0_solvents.itp" ]
        s_boxSize1 = s_boxSize #"6 6 6"
        s_boxSize2 = s_boxSize #"6.1 6.1 6.1"
        s_maxsol1 = s_maxsolW
        s_maxsol2 = s_maxsolO
    elif "charmm" in s_forceField.lower():
        solvent1   = "tip3p"
        solvent2   = [ solvents_dir / "charmm36" / "Octane" / "octane_box_npt.gro",  solvents_dir / "charmm36" / "Octane" / "octn.itp" ]
        s_boxSize1 = s_boxSize #"3 3 3"
        s_boxSize2 = s_boxSize #"3.1 3.1 3.1"
        s_maxsol1 = s_maxsolW
        s_maxsol2 = s_maxsolO





    # create the two systems it their own temporary folders
    s_molname = bricksFileSystem.get_filename_without_extension(s_pdbfile)
    pdb2molecule_in_solvent(s_pdbfile, s_molname + "_inW", solvent1, s_forceField, s_boxSize1, s_maxsol1, s_aditional_arguments)
    pdb2molecule_in_solvent(s_pdbfile, s_molname + "_inO", solvent2, s_forceField, s_boxSize2, s_maxsol2, s_aditional_arguments)

    #create folder that will contain the two systems, an moove them there
    bricksFileSystem.run_and_capture(f"mkdir {s_folderName}")
    bricksFileSystem.run_and_capture(f"rsync -av *_inW/* {s_folderName}")
    bricksFileSystem.run_and_capture(f"rsync -av *_inO/* {s_folderName}")
    bricksFileSystem.run_and_capture(f"rm -r *_inW")
    bricksFileSystem.run_and_capture(f"rm -r *_inO")




def create_tube_in_vacum(lipidsList, s_radius, s_thickness, s_box, s_outSysName, s_ff_location):
    """

    example:
    cl.create_tube_in_vacum([["DOPC","0.98","0.98","0.68"],["TO","0.02","0.02","0.68"]], s_radius= "18", s_thickness="2", s_box="50 50 50", s_outSysName="tube_in_vacum", s_ff_location="/home/bioinformatician/repositories/lipid_sorting/TS2CG-Setup-Pipeline/top/Martini3+NLs.LIB")    
    """

    # Create a folder to store the output system
    bricksFileSystem.run_and_capture(f"mkdir {s_outSysName}")
    s_outPathAndName = f"{s_outSysName}/{s_outSysName}"

    #copy the forcefiled to the system folder, and uptate the forcefield location
    bricksFileSystem.run_and_capture(f"cp -r \"{s_ff_location}\" \"{s_outSysName}\"")
    s_ff_location = os.path.basename(s_ff_location)
    s_ff_location = f"{s_outSysName}/{s_ff_location}"

    # Read the lipids list and format it so it can be pasted in the description file
    lipids_list_ready_to_paste = "\n".join(["     ".join(map(str, lipid)) for lipid in lipidsList])

    description_content = f"""
[Lipids List]
Domain 0
{lipids_list_ready_to_paste}
End

[Shape Data]
ShapeType Cylinder
Radius {s_radius}
Thickness {s_thickness}
Box {s_box}
End
    """

    # Save the description content to a temporary file
    with open("temporary_description.str", "w") as temp_file:
        temp_file.write(description_content)

    s_tool = "/home/bioinformatician/repositories/lipid_sorting/TS2CG-Setup-Pipeline/PCG"
    bricksFileSystem.run_and_capture(f"{s_tool} -str temporary_description.str -Bondlength 0.2 -LLIB {s_ff_location} -function analytical_shape -defout {s_outPathAndName}")

    # Delete the temporary file
    bricksFileSystem.delete("temporary_description.str")


def create_vesicle_in_vacum(lipidsList, s_radius, s_thickness, s_box, s_outSysName, s_ff_location):
    """

    example:
    cl.create_vesicle_in_vacum([["DOPC","0.98","0.98","0.68"],["TO","0.02","0.02","0.68"]], s_radius= "18", s_thickness="2", s_box="50 50 50", s_outSysName="vesicle_in_vacum", s_ff_location="/home/bioinformatician/repositories/lipid_sorting/TS2CG-Setup-Pipeline/top/Martini3+NLs.LIB")    
    """

    # Create a folder to store the output system
    bricksFileSystem.run_and_capture(f"mkdir {s_outSysName}")
    s_outPathAndName = f"{s_outSysName}/{s_outSysName}"

    #copy the forcefiled to the system folder, and uptate the forcefield location
    bricksFileSystem.run_and_capture(f"cp -r \"{s_ff_location}\" \"{s_outSysName}\"")
    s_ff_location = os.path.basename(s_ff_location)
    s_ff_location = f"{s_outSysName}/{s_ff_location}"

    # Read the lipids list and format it so it can be pasted in the description file
    lipids_list_ready_to_paste = "\n".join(["     ".join(map(str, lipid)) for lipid in lipidsList])

    description_content = f"""
[Lipids List]
Domain 0
{lipids_list_ready_to_paste}
End

[Shape Data]
ShapeType Sphere
Radius {s_radius}
Thickness {s_thickness}
Box {s_box}
End
    """

    # Save the description content to a temporary file
    with open("temporary_description.str", "w") as temp_file:
        temp_file.write(description_content)

    s_tool = "/home/bioinformatician/repositories/lipid_sorting/TS2CG-Setup-Pipeline/PCG"
    bricksFileSystem.run_and_capture(f"{s_tool} -str temporary_description.str -Bondlength 0.2 -LLIB {s_ff_location} -function analytical_shape -defout {s_outPathAndName}")

    # Delete the temporary file
    bricksFileSystem.delete("temporary_description.str")
