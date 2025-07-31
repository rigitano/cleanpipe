from cleanpipe import bricksPDB
from cleanpipe import bricksFileSystem
from cleanpipe import bricksTOP
from cleanpipe import bricksGROMACS

import subprocess
import re
import sys
import os
import functools

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


def pdb2molecule_in_solvent(s_pdbfile, s_outSytemName, s_solvent, s_forceField, s_boxSize):
    """
    s_pdbfile       : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
    s_outSytemName  : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    s_solvent       : choose a water model, for example as "tip3p", or a folder, for example "box_full_of_octn". The folder have to contain a system with a solvent box, in other words, it has to contain a octn_filledbox.gro and a octn.itp
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"

    example:
    cl.pdb2molecule_in_solvent("1LZ1.pdb", "1LZ1_in_water", "tip3p", "charmm36-jul2022", "3 3 3")
    cl.pdb2molecule_in_solvent("1LZ1.pdb", "1LZ1_in_octane", "box_full_of_octn", "charmm36-jul2022", "3 3 3")
    """

    # check if the filename inside s_pdbfile is valid
    bricksFileSystem.check_extention(s_pdbfile,['.pdb']) 

    # create gro, itps and top in a folder with the name of the system
    bricksGROMACS.pdb2system(s_pdbfile,s_outSytemName,s_forceField,s_boxSize)

    # add solvent to the system. I have 2 options here: tip3p or filled box
    bricksGROMACS.solvate_and_neutralize(s_outSytemName,s_solvent,s_forceField)

    # set the the name of the system in the top file 
    bricksTOP.setSystemName(f"{s_outSytemName}/{s_outSytemName}.top", f"{s_outSytemName} ; molecule from \"{s_pdbfile}\", inserted in solvent box of \"{s_solvent}\"" )



@ensure_original_directory
def pdb2molecule_in_two_solvents(s_pdbfile, s_folderName, s_outSytemName1, s_outSytemName2, s_solvent1, s_solvent2, s_forceField, s_boxSize):
    """
    s_pdbfile       : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
    s_folderName    : string with the name of the folder, for example "ala6Hdihr200-transfer". this folder will be created, and inside it, the two systems will be created inside it
    s_outSytemName1  : string with the name of the system, for example "ala6Hdih200_in_octane". a folder with that name will be created, and inside it, all the files, for example: ala6Hdihres200.gro and ala6Hdihres200.top
    s_outSytemName2  : string with the name of the system, for example "ala6Hdih200_in_water". a folder with that name will be created, and inside it, all the files, for example:ala6Hdihres200.gro and ala6Hdihres200.top
    s_solvent1       : for example "box_full_of_octn" folder with system with a solvent box, in other words, it has to contain a octn_filledbox.gro and a octn.itp
    s_solvent2       : for example as "tip3p"
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"

    example:
    cl.pdb2molecule_in_two_solvents("normal_peptide.pdb", "ala6Hdih200-transfer", "ala6Hdih200_in_octane", "ala6Hdih200_in_water", "box_full_of_octn", "tip3p", "charmm36-jul2022", "3 3 3")

    """

    # check if the filename inside s_pdbfile is valid
    bricksFileSystem.check_extention(s_pdbfile,['.pdb']) 

    
    #create output folder in parael with the pdb input. it will contain two folders. we will cd into each one of those do everithing there
    bricksFileSystem.create_folder(s_folderName)

    bricksFileSystem.run_and_capture(f"cp {s_pdbfile} {s_folderName}/temp.pdb")


    bricksFileSystem.run_and_capture(f"cp -r {s_forceField}.ff {s_folderName.rstrip('/')}/")#copy the forcefield to the new folder

    #I know in the example usage the first solvent is "box_full_of_octn", so I moove the folder. xxx this must be done differently
    bricksFileSystem.run_and_capture(f"cp -r {s_solvent1} {s_folderName.rstrip('/')}/")#copy the forcefield to the new folder

    try:
        bricksFileSystem.run_and_capture(f"cp -r toppar {s_folderName.rstrip('/')}/")#copy the forcefield to the new folder
    except:
        print("CLEANPIPE MESSAGE toppar folder not found, so it was not copied to the new folder")
    

    #go to the folder where the two systems will be created
    #original_directory = os.getcwd()#original folder is stored so I can go back to it at the very end of this function
    os.chdir(f"{s_folderName}")
    

    # CREATE FIRST SYSTEM (SOLVENT 1)
    bricksGROMACS.pdb2system("temp.pdb",s_outSytemName1,s_forceField,s_boxSize)   # create
    bricksGROMACS.solvate_and_neutralize(s_outSytemName1,s_solvent1,s_forceField) # solvate
    bricksTOP.setSystemName(f"{s_outSytemName1}/{s_outSytemName1}.top", f"{s_outSytemName1} ; molecule from \"{s_pdbfile}\", inserted in solvent box of \"{s_solvent1}\"" ) # set name

    # CREATE SECOND SYSTEM (SOLVENT 2)
    bricksGROMACS.pdb2system("temp.pdb",s_outSytemName2,s_forceField,s_boxSize)  # create
    bricksGROMACS.solvate_and_neutralize(s_outSytemName2,s_solvent2,s_forceField)# solvate
    bricksTOP.setSystemName(f"{s_outSytemName2}/{s_outSytemName2}.top", f"{s_outSytemName2} ; molecule from \"{s_pdbfile}\", inserted in solvent box of \"{s_solvent2}\"" ) # set name


    #remove files from the parante folder (the one that stores the two systems)
    bricksFileSystem.run_and_capture(f"rm temp.pdb")
    #I know in the example usage the first solvent is "box_full_of_octn", so I moove the folder. xxx this must be done differently
    bricksFileSystem.run_and_capture(f"rm -r {s_solvent1}")

    bricksFileSystem.run_and_capture(f"rm -r toppar")
    bricksFileSystem.run_and_capture(f"rm -r {s_forceField}.ff")

    #go back to the original folder
    #os.chdir(original_directory)


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
