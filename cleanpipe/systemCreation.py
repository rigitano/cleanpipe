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


import numpy as np
from scipy.spatial import cKDTree

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

@ensure_original_directory
def molecule2box_full_of_that(molecule, s_forceField, s_box_size, n_mol_max):
    """
    usage example: cl.molecule2box_full_of_that(["octn.gro", "octn.itp"], "charmm36-jul2022", "5 5 5", 1000)

    create a 5x5x5 box system filled with a lot of copies of the molecule
    as gromacs dont have such tool, its necessary to be creative, and use "insert-molecules"
    to fill the box with copyes of the molecule, then modify the top to reflect the new total
    the top file is also edited to change name of the system. and also the name of the molecule
    here the outputs are named box_full_of_blabla.gro and box_full_of_blabla.top, inside the
    folder box_full_of_blabla/
    """

    #if the variable 'molecule' contains a list, I will presume is a .gro and a .itp file names
    if isinstance(molecule, list) and len(molecule) == 2:

        #get the gro and itp as ABSOLUTE paths, so the os.chdir() done later cant break them
        s_gro = str(Path(molecule[0]).expanduser().resolve())
        s_itp = str(Path(molecule[1]).expanduser().resolve())

        #check if the filename inside s_gro and s_itp are valid
        bricksFileSystem.check_extention(s_gro, ['.gro'])
        bricksFileSystem.check_extention(s_itp, ['.itp'])

        #get the name of the molecule inside the itp
        ll_moleculetype = bricksTOP.parse_directive(s_itp, '[ moleculetype ]')
        s_extracted_mol_name = ll_moleculetype[0][0]

    else:
        raise ValueError("The 'molecule' parameter must be a list containing [gro_file, itp_file].")

    #the output folder must exist BEFORE anything is copied into it (this is why the molecule
    #is parsed first and the forcefield is copied only after)
    s_outName = f"box_full_of_{s_extracted_mol_name}"
    bricksFileSystem.run_and_capture(f'mkdir -p "{s_outName}"')  # -p so running it twice doesnt crash

    if str(s_forceField).lower() in ["charmm36-jul2022", "martini3001", "martini22"]:

        #if the user chose one the forcefields that I have stored myself in USEFUL_FORCEFIELDS
        module_path = Path(__file__).resolve().parent  # Where the cl module lives
        s_ffLocation = module_path / "USEFUL_FORCEFIELDS"

        if "charmm36-jul2022" in s_forceField.lower():
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/charmm36-jul2022.ff" "{s_outName}/"')  # in the case of charmm, the actual folder has a .ff in the end
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/toppar" "{s_outName}/"')  # and this extra file must alse come
            ff_inclusion_text = "charmm36-jul2022.ff/forcefield.itp"

        elif "martini3001" in s_forceField.lower():
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/martini3001" "{s_outName}/"')
            ff_inclusion_text = "martini3001/martini_v3.0.0.itp"

        elif "martini22" in s_forceField.lower():
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/martini22" "{s_outName}/"')
            ff_inclusion_text = "martini22/martini_v2.2.itp"
    else:

        #if the user gave the folder of the forcefield
        s_forceField = str(Path(s_forceField).expanduser().resolve())  #resolve all ../ ~/ ../../ ./ to absolute path
        s_ff_foldername = Path(s_forceField).name  # the folder name AS IT IS, a .ff folder must keep its .ff
        bricksFileSystem.run_and_capture(f'cp -r "{s_forceField}" "{s_outName}/"')
        ff_inclusion_text = f"{s_ff_foldername}/forcefield.itp"

    ############### manualy create top ###############

    s_itp_name = Path(s_itp).name  # the top will sit next to the copyed itp, so only the filename goes in the #include
    topology_text = f"""
#include "{ff_inclusion_text}"
#include "{s_itp_name}"

[ system ]
box full of {s_extracted_mol_name}

[ molecules ]
{s_extracted_mol_name}   1
"""
    with open(f"{s_outName}/{s_outName}.top", "w") as f:  # written directly inside the folder, no need to mv it later
        f.write(topology_text)

    ###################################################

    # bring the original gro and itp to the system folder
    bricksFileSystem.run_and_capture(f'cp "{s_gro}" "{s_outName}/"')
    bricksFileSystem.run_and_capture(f'cp "{s_itp}" "{s_outName}/"')

    #cd into the output folder we created and do everithing there
    os.chdir(s_outName)

    #manipulate the GRO file to create a box and fill it with copyes of the molecule
    #(from here on the names are relative to the folder we are already inside of)
    captured_output = bricksFileSystem.run_and_capture(f'gmx insert-molecules -ci "{Path(s_gro).name}" -nmol {str(n_mol_max)} -rot xyz -box {s_box_size} -o "{s_outName}.gro"')
    print(f"\nCLEANPIPE MESSAGE\ngro file written: \n {s_outName}/{s_outName}.gro")

    #get the number of added molecules.
    match = re.search(r'Added\s+(\d+)\s+molecules', captured_output)
    if match is None:  # if gmx failed or printed something else, dont crash with a cryptic NoneType error
        raise RuntimeError(f"could not read how many molecules 'gmx insert-molecules' added:\n{captured_output}")
    added_molecules = int(match.group(1))

    #update the TOP file with the new total the molecule
    bricksTOP.update_molecule_quantity(f"{s_outName}.top", s_extracted_mol_name, added_molecules)
    

@ensure_original_directory
def molecule2molecule_in_solvent(molecule, s_outSytemName, solvent, s_forceField, s_boxSize, s_maxsol=0, s_aditional_arguments=''):
    """
    molecule       : string with the pdb name. for example "insulin.pdb", 
                      or list with gro and itp names. for example ["insulin.gro", "insulin.itp"]
                        this will be the main molecule in the system.
    s_outSytemName  : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    solvent         : choose a water model, for example as "tip3p", or a list containg a gro of a box of solvents and its itp, for example ["../solvents/box_full_of_octn.gro","octn.itp], bot the file name should be set in reference to the system top
    s_forceField    : place where the ff is, relative to the current location, so to be copyed to the system folder. examples:
                                                                                                                            '../alanine12/v15-truss/charmm36-jul2022.ff'
                                                                                                                            '~/ff/martini3001'
    s_boxSize       : string with x y z sizes, for example "3 3 3"
    s_maxsol        : the maximum number of solvent molecules that will be added. this is optional here, as the 0 value mean the parameter wont be considered by gromacs

    s_aditional_arguments : add extra arguments in martinize or pdb2gmx. for example, in martinize: "-water-bias -water-bias-eps E:-0.5 H:-1.0 -ss HHHHHHHHHHHH"

    examples:

    cl.molecule2molecule_in_solvent("g12H.pdb", "g12HW", "tip3p", 'charmm36-jul2022', "6 6 6", "6943")
    cl.molecule2molecule_in_solvent("g12H.pdb", "g12HO", ["~/repos/cleanpipe/cleanpipe/USEFUL_SOLVENT_BOXES/charmm36/Wet_Octanol/wet_octanol.gro","~/repos/cleanpipe/cleanpipe/USEFUL_SOLVENT_BOXES/charmm36/Octane/octn.itp"], 'charmm36-jul2022', "6.1 6.1 6.1", "712")

    cl.molecule2molecule_in_solvent("g12H.pdb", "g12HW_cg", ["~/solvents/martini/Water-pure/water.gro","~/ff/martini3001/martini_v3.0.0_solvents_v1.itp"], '~/ff/martini3001', "6 6 6", "1736", "-ss HHHHHHHHHHHH")
    cl.molecule2molecule_in_solvent("g12H.pdb", "g12HO_cg", ["~/solvents/martini/Octane/OCT_PRO1.gro", "~/ff/martini3001/martini_v3.0.0_solvents_v1.itp"], '~/ff/martini3001', "6.1 6.1 6.1", "712", "-ss HHHHHHHHHHHH")

    """
    print("CLEANPIPE called function molecule2molecule_in_solvent")


    #if the variable 'molecule' contains a list, I will presume is a .gro and a .top file names
    if isinstance(molecule, list) and len(molecule) == 2:
        #get the gro and top
        s_gro = molecule[0]
        s_origin = molecule[1]

        #check if the filename inside s_gro and s_itp are valid
        bricksFileSystem.check_extention(s_gro,['.gro']) 
        bricksFileSystem.check_extention(s_origin,['.itp']) 


    else: # I will presume its a pdb filename

        s_origin = molecule

        #check if the filename is valid
        bricksFileSystem.check_extention(s_origin,['.pdb']) 







    ########################### create gro, itps and top in a folder with the name of the system ###########################
    bricksGROMACS.molecule2system(molecule, s_outSytemName, s_forceField, s_boxSize, s_aditional_arguments)
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
    bricksTOP.setSystemName(f"{s_outSytemName}/{s_outSytemName}.top", f"{s_outSytemName} ; molecule from \"{s_origin}\", inserted in solvent from \"{s_solvent_text}\"" )


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
def molecule2molecule_in_water_and_oil(molecule, s_oil_choice, s_folderName, s_forceField, s_boxSize ,s_maxsolW, s_maxsolO, s_aditional_arguments=''):
    """
    this function is perfect for canclulations of free energies of transfer. there will be two systems in the same folder. 
    but for DG transfer, this actualy makes things easyer. you can run two runREALISTIC scrits, and then two runFEPoff scripts. and all the data you need will be there

    molecule        : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
                       or list with gro and itp names. for example ["insulin.gro", "insulin.itp"] 
    s_oil_choice    : for now, you have to chose among the possible options, because inside this function, I have to hardcode the adress of the box of solvent:
                            octane
                            wet_octanol (just for charmm36)
                      
    
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
    cl.molecule2molecule_in_water_and_oil("normal_peptide.pdb", "octane", "ala6Hdih200-transfer", "charmm36-jul2022", "6 6 6", "6943", "712" )
    cl.molecule2molecule_in_water_and_oil("g12H.pdb",           "octane", "g12Hbias1-transfer",   '~/ff/martini3001', "6 6 6", "1736", "712", s_aditional_arguments='-water-bias -water-bias-eps E:-0.5 H:-1.0 -ss HHHHHHHHHHHH')


    """
    print("CLEANPIPE called function molecule2molecule_in_water_and_oil")




    #if the variable 'molecule' contains a list, I will presume is a .gro and a .top file names
    if isinstance(molecule, list) and len(molecule) == 2:
        #get the gro and top
        s_gro = molecule[0]
        s_origin = molecule[1]

        #check if the filename inside s_gro and s_itp are valid
        bricksFileSystem.check_extention(s_gro,['.gro']) 
        bricksFileSystem.check_extention(s_origin,['.itp']) 


    else: # I will presume its a pdb filename

        s_origin = molecule

        #check if the filename is valid
        bricksFileSystem.check_extention(s_origin,['.pdb']) 








    module_path = Path(__file__).resolve().parent # Where the cl module lives
    ff_dir       = module_path / "USEFUL_FORCEFIELDS" 
    solvents_dir = module_path / "USEFUL_SOLVENT_BOXES" 

    #get the water
    if "martini3001" in s_forceField.lower():
        solvent1  = [ solvents_dir / "martini3001" / "Water-pure" / "water.gro",     ff_dir / "martini3001" / "martini_v3.0.0_solvents_v1.itp" ]    # I think water is already in the ff itp, so no need for this solvents file       

    elif "martini22" in s_forceField.lower():
        solvent1  = [ solvents_dir / "martini22" / "Water-pure" / "water.gro",     ff_dir / "martini22" / "martini_v2.0_solvents.itp" ] # I think water is already in the ff itp, so no need for this solvents file  

    elif "charmm" in s_forceField.lower():
        solvent1   = "tip3p"

    #get the oil
    if s_oil_choice == "octane":
        if "martini3001" in s_forceField.lower():
            solvent2  = [ solvents_dir / "martini3001" / "Octane"     / "OCT_PRO1.gro",  ff_dir / "martini3001" / "martini_v3.0.0_solvents_v1.itp" ]

        elif "martini22" in s_forceField.lower():
            solvent2  = [ solvents_dir / "martini22" / "Octane"     / "OCT_PRO1.gro",  ff_dir / "martini22" / "martini_v2.0_solvents.itp" ]

        elif "charmm" in s_forceField.lower():
            solvent2   = [ solvents_dir / "charmm36" / "Octane" / "octane_box_npt.gro",  solvents_dir / "charmm36" / "Octane" / "octn.itp" ]

    elif s_oil_choice == "wet_octanol":
        if "martini3001" in s_forceField.lower():
            solvent2  = [ solvents_dir / "martini3001" / "xxxxx"     / "xxxxx",  ff_dir / "martini3001" / "xxxxx" ]

        elif "martini22" in s_forceField.lower():
            solvent2  = [ solvents_dir / "martini22" / "xxxxx"     / "xxxxxx",  ff_dir / "martini22" / "xxxxxxx" ]

        elif "charmm" in s_forceField.lower():
            solvent2   = [ solvents_dir / "charmm36" / "Wet_Octanol" / "prod.gro",  solvents_dir / "charmm36" / "Wet_Octanol" / "OCTO.itp",  solvents_dir / "charmm36" / "Wet_Octanol" / "tip3p.itp" ]




    s_boxSize1 = s_boxSize #for example "3 3 3"
    s_boxSize2 = s_boxSize #for example "3 3 3"
    s_maxsol1 = s_maxsolW
    s_maxsol2 = s_maxsolO


    # create the two systems it their own temporary folders
    s_molname = bricksFileSystem.get_filename_without_extension(s_origin)
    molecule2molecule_in_solvent(molecule, s_molname + "_inW", solvent1, s_forceField, s_boxSize1, s_maxsol1, s_aditional_arguments)
    molecule2molecule_in_solvent(molecule, s_molname + "_inO", solvent2, s_forceField, s_boxSize2, s_maxsol2, s_aditional_arguments)

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



@ensure_original_directory
def build_membrane(
    out_system_name,
    lipids_by_domain,
    shape_type,
    shape_params,
    proteins={},
    floating_domains=None,
    thickness=3.8,
    density=3,
    bondlength=0.2,
    llib='Martini3.LIB',
    
    default_exclusion_domain=0
):
    """
    Automated TS2CG workflow:
    1. Writes an input .str file including lipids and proteins.
    2. Generates a point folder using PCG with -WPointDir.
    3. Sequentially inserts proteins using INU.
    4. Optionally assigns lipid domains around proteins (DAI).
    5. Optionally excludes lipids from pore proteins (DAI).
    6. Optionally creates floating lipid domains (DAI with manual points).
    7. Builds the final .gro/.top system with PCG.

    Parameters
    ----------


    out_system_name : str
        Name of the folder that will be created, and also the .gro and .top that will appear there).
    lipids_by_domain : dict
        Mapping of domain IDs to lipid definitions.
        Format: {domain: [(lipid_name, ratio_up, ratio_down, APL), ...]}.
    shape_type : str
        One of {"Sphere", "Cylinder", "Flat", "1D_PBC_Fourier"}.
    shape_params : dict
        Shape-specific parameters (e.g. {"Box": (60,60,60), "Radius": 25}).
    proteins : dict
        Protein definitions.
        Example:
        {
          "protein_1": {
              "gro": "../structures/protein_1.gro",
              "type_id": 1,
              "count": 2,
              "radius": 5,
              "pore": False,
              "domain_assignment": {"domain_id": 1, "radius": 5}
          },
          "protein_2": {
              "gro": "../structures/protein_2.gro",
              "type_id": 2,
              "count": 3,
              "radius": 5,
              "pore": True,
              "pore_radius": 6
          }
        }
    floating_domains : list of dict, optional
        Floating circular lipid domains not tied to proteins.
        Example:
        [{"domain_id": 1, "radius": 8, "points": [200, 450]}]
        - domain_id must exist in [Lipids List].
        - points are vertex indices from the point file.
    thickness : float
        Bilayer thickness (nm).
    density : float or tuple
        Point density for analytical shapes.

    bondlength : float
        Bond length for PCG.

    default_exclusion_domain : int
        Lipid domain ID to assign when making pore exclusions.
        Must exist in [Lipids List].
    llib : str
        Path to Martini3.LIB file. example "../Martini3.LIB". 
        be carefull, the Martini3.LIB from USEFUL_LISTS will be copied to the current folder, so if you want to use a different one, you have to provide a different path to it.
    """


    Path(out_system_name).mkdir(parents=True, exist_ok=True)
    os.chdir(out_system_name)


    out_specs_filename = "system_specification.str"

    # ------------------------------------------------------------------
    # Helper function: run command with error checking
    # ------------------------------------------------------------------
    def run_cmd(cmd, desc):
        print(f"\n--- {desc} ---")
        print("Command:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            sys.exit(f"ERROR during {desc}. Command failed:\n{e}\n"
                     f"Check TS2CG output above for details.")
        print(f"Step completed successfully: {desc}")

    # ------------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------------
    # 1. Check protein files exist
    for pname, pinfo in proteins.items():
        if not os.path.exists(pinfo["gro"]):
            sys.exit(f"ERROR: Protein file {pinfo['gro']} for {pname} not found.")

    # 2. Check lipid library exists after bringing it from USEFUL_LISTS
    module_path = Path(__file__).resolve().parent
    s_llLocation = module_path / "USEFUL_LISTS"
    bricksFileSystem.run_and_capture(f'cp "{s_llLocation}/Martini3.LIB" .')

    if not os.path.exists(llib):
        sys.exit(f"ERROR: Lipid library file {llib} not found. "
                 f"Provide the correct path to Martini3.LIB.")

    # 3. Check lipid ratios are reasonable (sum ≈ 1)
    for domain, lipids in lipids_by_domain.items():
        total_up = sum(l[1] for l in lipids)
        total_down = sum(l[2] for l in lipids)
        if abs(total_up - 1.0) > 0.2 or abs(total_down - 1.0) > 0.2:
            print(f"WARNING: Domain {domain} lipid ratios may be inconsistent "
                  f"(sum_up={total_up:.2f}, sum_down={total_down:.2f}).")

    # 4. Check default exclusion domain exists
    if default_exclusion_domain not in lipids_by_domain:
        sys.exit(f"ERROR: default_exclusion_domain {default_exclusion_domain} "
                 f"is not defined in [Lipids List].")

    # ------------------------------------------------------------------
    # Step 1: Write .str file
    # ------------------------------------------------------------------
    print("\nWriting input .str file:", out_specs_filename)
    with open(out_specs_filename, "w") as f:
        f.write("; Auto-generated TS2CG input\n")

        # Include protein coordinate files
        for pname, pinfo in proteins.items():
            f.write(f"include {pinfo['gro']}\n")
        f.write("\n")

        # Lipid domains
        f.write("[Lipids List]\n")
        for domain, lipids in lipids_by_domain.items():
            f.write(f"Domain   {domain}\n")
            for name, up, down, apl in lipids:
                f.write(f"{name:<8} {up:<6} {down:<6} {apl:<6}\n")
            f.write("End\n\n")

        # Protein list
        if proteins:
            f.write("[Protein List]\n")
            for pname, pinfo in proteins.items():
                f.write(f"{pname:<12} {pinfo['type_id']} 0.01 0 0 -2.0\n")
            f.write("End Protein\n\n")

        # Shape
        f.write("[Shape Data]\n")
        f.write(f"ShapeType   {shape_type}\n")
        if "Box" in shape_params:
            f.write("Box   " + "   ".join(map(str, shape_params["Box"])) + "\n")
        if density is not None:
            if isinstance(density, (tuple, list)):
                f.write("Density   " + "   ".join(map(str, density)) + "\n")
            else:
                f.write(f"Density   {density}\n")
        f.write(f"Thickness   {thickness}\n")

        if shape_type == "Sphere":
            if "WallDensity" in shape_params:
                f.write("WallDensity   " + "   ".join(map(str, shape_params["WallDensity"])) + "\n")
            if "DL" in shape_params:
                f.write(f"DL   {shape_params['DL']}\n")
            f.write(f"Radius   {shape_params['Radius']}\n")
        elif shape_type == "Cylinder":
            f.write(f"Radius   {shape_params['Radius']}\n")
        elif shape_type == "Flat":
            if "WallRange" in shape_params:
                f.write("WallRange   " + "   ".join(map(str, shape_params["WallRange"])) + "\n")
        elif shape_type == "1D_PBC_Fourier":
            if "WallRange" in shape_params:
                f.write("WallRange   " + "   ".join(map(str, shape_params["WallRange"])) + "\n")
            if "Modes" in shape_params:
                for mode in shape_params["Modes"]:
                    f.write("Mode   " + "   ".join(map(str, mode)) + "\n")
        f.write("End\n")
    print("Input .str file written successfully.")

    # ------------------------------------------------------------------
    # Step 2: Generate initial point folder
    # ------------------------------------------------------------------
    if os.path.exists("point"):
        shutil.rmtree("point")
    run_cmd([
        "TS2CG", "PCG",
        "-str", out_specs_filename,
        "-Bondlength", str(bondlength),
        "-LLIB", llib,
        "-function", "analytical_shape",
        "-defout", out_system_name,
        "-WPointDir"
    ], "Generate initial point folder")

    # ------------------------------------------------------------------
    # Step 3: Insert proteins with INU, apply domain assignments and exclusions
    # ------------------------------------------------------------------
    point_dir = "point"
    for i, (pname, pinfo) in enumerate(proteins.items(), 1):
        outdir = f"point_new{i}"

        # Insert proteins
        run_cmd([
            "TS2CG", "INU",
            "--point-dir", point_dir,
            "--protein-type", str(pinfo["type_id"]),
            "--radius", str(pinfo["radius"]),
            "--num-proteins", str(pinfo["count"]),
            "--output-dir", outdir
        ], f"Insert {pinfo['count']} copies of {pname}")
        point_dir = outdir

        # Optional: assign lipid domain around protein
        if "domain_assignment" in pinfo:
            da = pinfo["domain_assignment"]
            domain_out = f"{outdir}_dom"
            run_cmd([
                "TS2CG", "DAI",
                "--point-dir", outdir,
                "--radius", str(da["radius"]),
                "--domain-id", str(da["domain_id"]),
                "--protein-type", str(pinfo["type_id"]),
                "--output-dir", domain_out
            ], f"Assign domain {da['domain_id']} around {pname}")
            point_dir = domain_out

        # Optional: pore exclusion
        if pinfo.get("pore", False):
            excl_outdir = f"{point_dir}_excl"
            pore_r = pinfo.get("pore_radius", pinfo["radius"])
            run_cmd([
                "TS2CG", "DAI",
                "--point-dir", point_dir,
                "--radius", str(pore_r),
                "--domain-id", str(default_exclusion_domain),
                "--protein-type", str(pinfo["type_id"]),
                "--output-dir", excl_outdir
            ], f"Add lipid exclusion for pore protein {pname}")
            point_dir = excl_outdir

    # ------------------------------------------------------------------
    # Step 3b: Floating domains (not tied to proteins)
    # ------------------------------------------------------------------
    if floating_domains:
        for j, dom in enumerate(floating_domains, 1):
            outdir = f"{point_dir}_float{j}"
            run_cmd([
                "TS2CG", "DAI",
                "--point-dir", point_dir,
                "--radius", str(dom["radius"]),
                "--domain-id", str(dom["domain_id"]),
                "--manual-points", ",".join(map(str, dom["points"])),
                "--output-dir", outdir
            ], f"Assign floating domain {dom['domain_id']} at points {dom['points']}")
            point_dir = outdir

    # ------------------------------------------------------------------
    # Step 4: Build final system
    # ------------------------------------------------------------------
    run_cmd([
        "TS2CG", "PCG",
        "-str", out_specs_filename,
        "-Bondlength", str(bondlength),
        "-LLIB", llib,
        "-dts", point_dir,
        "-defout", out_system_name
    ], "Build final system")

    print(f"\nFinal system successfully built: {out_system_name}.gro and {out_system_name}.top")


    # Bring Martini force field to the folder
    module_path = Path(__file__).resolve().parent
    s_ffLocation = module_path / "USEFUL_FORCEFIELDS"
    bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/martini3001" .')



    # Add Martini includes at the beginning of the topology
    top_file = Path(f"{out_system_name}.top")

    martini_includes = (
        '#include "martini3001/martini_v3.0.0.itp"\n'
        '#include "martini3001/martini_v3.0.0_ffbonded_v2.itp"\n'
        '#include "martini3001/martini_v3.0.0_phospholipids_v1.itp"\n'
        '#include "martini3001/martini_v3.0_sterols_v1.0.itp"\n'
        '#include "martini3001/martini_v3.0.0_solvents_v1.itp"\n'
        '#include "martini3001/martini_v3.0.0_ions_v1.itp"\n'
        '#include "martini3001/martini_v3.0.0_sugars_v1.itp"\n'
        '\n'
    )

    top_file.write_text(
        martini_includes + top_file.read_text()
    )

    files_to_copy = [
        "runREALISTIC.sh",
        "runBENCHMARK-rome.sh",
    ]
    module_path = Path(__file__).resolve().parent # Where the cl module lives
    source_dir = module_path / "bash" # Folder containing the source files
    dest_dir = Path(out_system_name) # Destination folder
    for filename in files_to_copy:
        bricksFileSystem.run_and_capture(f'cp -r "{source_dir}/{filename}" .')





def slab_in_water(gro_in,
                        top_in,
                        layer_thickness,          # nm, per layer
                        out_dir,
                        s_forceField,
                        solvent_box="spc216.gro",  # pre-equilibrated water box
                        min_dist=0.22,             # nm, water<->slab clash cutoff
                        gmx="gmx",
                        keep_workdir=False,
                        verbose=True):

    """


    Build a  [ water | your slab | water ]  sandwich along z, for interfacial
    tension calculations, WITHOUT ever inserting water inside your original system.

    Side view (z is vertical), Lz = original box height, t = layer thickness:

        2t + Lz  +---------------------+  <- new box top  (PBC partner of the bottom)
                |    water layer B    |   t
        t + Lz  +---------------------+
                |                     |
                |   your original     |   Lz   (atoms untouched, only shifted up)
                |   .gro system       |
            t  +---------------------+
                |    water layer A    |   t
            0  +---------------------+  <- new box bottom


    Why not "make the box taller and run gmx solvate on the whole thing"?
    Because solvate would also push water into every cavity of your slab.
    Here the water is generated on its own and merely *stacked* around the slab.


    The trick used to make the water (important, read this):
    --------------------------------------------------------
    * ONE pre-equilibrated water box of size (Lx, Ly, 2t) is generated with
    `gmx solvate -cs spc216.gro -box Lx Ly 2t`.
    * It is cut in half at z = t. The lower half becomes layer A, the upper half
    becomes layer B (shifted up by Lz).
    -> The two faces that end up meeting across the *outer* periodic boundary are
        exactly the two faces that were already periodic partners inside the
        original water box. So that junction is perfectly packed, and the two
        layers are still two genuinely different water configurations (not copies).
    * The two *inner* faces touch your slab. Waters that bump into it are deleted
    whole-molecule, which is exactly what `gmx solvate` does when it solvates a
    protein.

    Outputs go into a brand-new folder. The inputs are never touched.

    Dependencies: numpy, scipy, and your own `bricksFileSystem` helper.
    """


    """
    Put a water layer of `layer_thickness` nm above AND below the system in
    `gro_in`, keeping x and y exactly as they are.

    Parameters
    ----------
    gro_in, top_in   : paths to the original system (never modified)
    layer_thickness  : thickness in nm of EACH water layer
    out_dir          : new folder to create; everything is written there
    out_name         : basename of the outputs -> out_name.gro / out_name.top
    solvent_box      : solvent coordinates for `gmx solvate -cs`
                       (spc216.gro is fine for any 3-point model: SPC, SPC/E, TIP3P)

    min_dist         : a whole water molecule is deleted if any of its atoms is
                       closer than this to any slab atom (0.22 nm ~ what
                       `gmx solvate` uses by default for C/N/O)
    gmx              : name/path of the GROMACS binary
    keep_workdir     : keep the intermediate files for debugging
    verbose          : print a short report

    Returns
    -------
    dict with 'gro', 'top', 'n_water', 'box'


    example:
    cl.slab_in_water(
        gro_in="slab.gro",
        top_in="slab.top",
        s_forceField="charmm36-jul2022",
        layer_thickness=3.0,              # nm of water on each side
        out_dir="slab_test",
    )

    """




    def _keep_non_clashing(wat_xyz, wat_mol, slab_xyz, box, min_dist):
        """Boolean mask over water atoms: False for every atom of a water molecule
        that has at least one atom closer than `min_dist` to a slab atom.

        The KD-tree is built with `boxsize=box`, i.e. it measures distances through
        the periodic boundaries in x, y and z - which is what the real simulation
        will do too.
        """
        # KD-tree with PBC needs every coordinate inside [0, L)
        wrap = lambda p: np.where((p % box) >= box, 0.0, (p % box))

        tree = cKDTree(wrap(slab_xyz), boxsize=box)
        dist, _ = tree.query(wrap(wat_xyz), k=1, distance_upper_bound=min_dist)
        clashing_atom = np.isfinite(dist)        # inf = nothing within min_dist

        bad_molecules = np.unique(wat_mol[clashing_atom])
        return ~np.isin(wat_mol, bad_molecules)

    # =====================================================================
    # 1. Minimal .gro reader / writer
    #    .gro is a FIXED-COLUMN text format, all lengths in nm:
    #       cols  0:5   residue number
    #       cols  5:10  residue name
    #       cols 10:15  atom name
    #       cols 15:20  atom number
    #       then 3 (or 6, if velocities) equal-width float fields
    #    Last line = box vectors.
    # =====================================================================

    def _coord_field_width(atom_line):
        """Figure out how wide the coordinate columns are (8 by default, but files
        written with extra precision use wider fields)."""
        body = atom_line[20:].rstrip()
        for n_fields in (3, 6):                      # xyz, or xyz + velocities
            if len(body) % n_fields == 0:
                w = len(body) // n_fields
                if 7 <= w <= 15:                     # sane column width
                    return w
        return 8                                     # GROMACS default


    def _read_gro(path):
        """Read a .gro file. Velocities are ignored (they are meaningless after we
        rebuild the system anyway; minimisation/equilibration will regenerate them)."""
        with open(path) as fh:
            lines = fh.read().splitlines()

        title    = lines[0]
        n_atoms  = int(lines[1].strip())
        body     = lines[2:2 + n_atoms]
        box_line = lines[2 + n_atoms]

        w = _coord_field_width(body[0])

        resid    = np.empty(n_atoms, dtype=int)
        resname  = []
        atomname = []
        xyz      = np.empty((n_atoms, 3), dtype=float)

        for i, line in enumerate(body):
            resid[i] = int(line[0:5])
            resname.append(line[5:10].strip())
            atomname.append(line[10:15].strip())
            # atom number (cols 15:20) is dropped: we renumber everything on output
            c = line[20:]
            xyz[i] = [float(c[k * w:(k + 1) * w]) for k in range(3)]

        box = np.array([float(v) for v in box_line.split()], dtype=float)

        return dict(title=title, resid=resid, resname=resname,
                    atomname=atomname, xyz=xyz, box=box)


    def _write_gro(path, title, resid, resname, atomname, xyz, box_xyz):
        """Write a rectangular-box .gro file with the standard GROMACS formatting."""
        n = len(xyz)
        with open(path, "w") as fh:
            fh.write(title.strip() + "\n")
            fh.write("%d\n" % n)
            for i in range(n):
                fh.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (
                    resid[i] % 100000,            # .gro fields wrap at 99999
                    resname[i][:5],
                    atomname[i][:5],
                    (i + 1) % 100000,
                    xyz[i, 0], xyz[i, 1], xyz[i, 2]))
            fh.write("%10.5f%10.5f%10.5f\n" % tuple(box_xyz))


    def _molecule_ids(resid):
        """Give every atom the index of the molecule it belongs to.
        A new molecule starts whenever the residue number changes (true for water:
        GROMACS requires one residue per solvent molecule)."""
        is_new = np.ones(len(resid), dtype=bool)
        is_new[1:] = resid[1:] != resid[:-1]
        return np.cumsum(is_new) - 1


    # =====================================================================
    # 2. Topology helper: add the new water molecules to the .top
    # =====================================================================

    def _patch_topology(top_path, water_resname, n_water, water_itp_include=None):
        """Append '<water_resname>  <n_water>' to [ molecules ], and optionally add
        the #include line for the water model.

        NOTE: [ molecules ] must list molecules in the SAME ORDER as the atoms in
        the .gro. We write the solute first and all the water last, so appending at
        the end of the file is correct.
        """
        with open(top_path) as fh:
            lines = fh.read().splitlines()


        

        # -- optionally insert the water model #include right after the force field
        if water_itp_include:
            already_there = any(water_itp_include in l for l in lines)
            if not already_there:
                # put it just after the force field include, else after the 1st include
                pos = None
                for i, l in enumerate(lines):
                    if l.strip().startswith("#include") and "forcefield.itp" in l:
                        pos = i + 1
                        break
                if pos is None:
                    for i, l in enumerate(lines):
                        if l.strip().startswith("#include"):
                            pos = i + 1
                            break
                if pos is None:
                    pos = 0
                lines.insert(pos, '#include "%s"' % water_itp_include)

        # -- append the solvent count at the very end (= end of [ molecules ])
        lines.append("%-15s %d" % (water_resname, n_water))

        with open(top_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")


    # the main code

    out_name=out_dir          # basename of the outputs -> out_name.gro / out_name.top


    # ---------------------------------------------------------------
    # 3.1 Read the original system and sanity-check it
    # ---------------------------------------------------------------
    gro_in = os.path.abspath(gro_in)
    top_in = os.path.abspath(top_in)
    out_dir = os.path.abspath(out_dir)

    slab = _read_gro(gro_in)
    box = slab["box"]

    # Only rectangular boxes make sense for a flat slab geometry.
    if len(box) > 3 and np.any(np.abs(box[3:]) > 1e-6):
        raise ValueError("The input box is triclinic. A slab needs a rectangular box "
                         "(use `gmx editconf -bt cubic` / -box first).")
    Lx, Ly, Lz = box[0], box[1], box[2]

    if abs(Lx - Ly) > 1e-3:
        print("[warning] x (%.3f) and y (%.3f) differ - you said 'square', "
              "double-check this is what you want." % (Lx, Ly))
    if layer_thickness <= 0:
        raise ValueError("layer_thickness must be > 0 nm")
    if layer_thickness < 1.0:
        print("[warning] a %.2f nm layer is thin: the two interfaces will feel "
              "each other through PBC. 2-3 nm is a safer minimum." % layer_thickness)

    t = float(layer_thickness)
    Lz_new = Lz + 2.0 * t
    new_box = np.array([Lx, Ly, Lz_new])

    # ---------------------------------------------------------------
    # 3.2 Create the output folder + a scratch subfolder, copy the topology
    # ---------------------------------------------------------------
    work_dir = os.path.join(out_dir, "_work")
    bricksFileSystem.run_and_capture(f'mkdir -p "{work_dir}"')

    top_out = os.path.join(out_dir, out_name + ".top")
    gro_out = os.path.join(out_dir, out_name + ".gro")

    # copy the .top ... and any .itp sitting next to it (they are usually
    # #included with a relative path, so they must travel with the .top)
    bricksFileSystem.run_and_capture(f'cp "{top_in}" "{top_out}"')
    src_dir = os.path.dirname(top_in)
    for f in sorted(os.listdir(src_dir)):
        if f.endswith(".itp"):
            bricksFileSystem.run_and_capture(
                f'cp "{os.path.join(src_dir, f)}" "{os.path.join(out_dir, f)}"')

    # ---------------------------------------------------------------
    # 3.3 Generate ONE water box of size (Lx, Ly, 2t)  -- pure solvent,
    #     no -cp, so nothing of yours is involved.
    # ---------------------------------------------------------------
    water_gro = os.path.join(work_dir, "water_box.gro")
    bricksFileSystem.run_and_capture(
        f'{gmx} solvate -cs {solvent_box} '
        f'-box {Lx:.5f} {Ly:.5f} {2.0 * t:.5f} '
        f'-o "{water_gro}"')

    water = _read_gro(water_gro)
    w_xyz = water["xyz"]
    w_mol = _molecule_ids(water["resid"])           # which atoms form one molecule
    water_resname = water["resname"][0]             # normally "SOL"

    # ---------------------------------------------------------------
    # 3.4 Cut the water box in half at z = t and move the two halves into place
    #     A molecule goes with its centre, so molecules are never chopped.
    # ---------------------------------------------------------------
    n_mol = w_mol.max() + 1
    # mean z of each molecule (bincount = fast "group by molecule and average")
    mol_z = np.bincount(w_mol, weights=w_xyz[:, 2]) / np.bincount(w_mol)
    mol_is_lower = mol_z < t                        # lower half -> layer A
    atom_is_lower = mol_is_lower[w_mol]

    lower_xyz = w_xyz[atom_is_lower].copy()                 # stays at z in [0, t]
    upper_xyz = w_xyz[~atom_is_lower].copy()
    upper_xyz[:, 2] += Lz                                   # pushed above the slab

    lower_mol = w_mol[atom_is_lower]
    upper_mol = w_mol[~atom_is_lower]

    # ---------------------------------------------------------------
    # 3.5 Lift the original system by t, so it sits between the two layers
    # ---------------------------------------------------------------
    slab_xyz = slab["xyz"].copy()
    slab_xyz[:, 2] += t

    # ---------------------------------------------------------------
    # 3.6 Delete water molecules that clash with the slab
    #     (whole molecules only, otherwise the topology breaks)
    # ---------------------------------------------------------------
    wat_xyz = np.vstack([lower_xyz, upper_xyz])
    # relabel molecules 0..N-1 over the two concatenated halves
    wat_mol = np.concatenate([lower_mol, upper_mol + n_mol])
    _, wat_mol = np.unique(wat_mol, return_inverse=True)

    keep_atom = _keep_non_clashing(wat_xyz, wat_mol, slab_xyz, new_box, min_dist)

    wat_xyz = wat_xyz[keep_atom]
    wat_mol = wat_mol[keep_atom]
    _, wat_mol = np.unique(wat_mol, return_inverse=True)     # renumber again
    n_water = int(wat_mol.max()) + 1 if len(wat_mol) else 0

    # water atom/residue names, in the same kept order
    w_resname_all = np.array(water["resname"])
    w_atomname_all = np.array(water["atomname"])
    order = np.concatenate([np.where(atom_is_lower)[0], np.where(~atom_is_lower)[0]])
    wat_resname = list(w_resname_all[order][keep_atom])
    wat_atomname = list(w_atomname_all[order][keep_atom])

    # ---------------------------------------------------------------
    # 3.7 Write the merged .gro:  SLAB FIRST, then all the water.
    #     (order must match [ molecules ] in the .top)
    # ---------------------------------------------------------------
    first_water_resid = int(slab["resid"].max()) + 1
    all_resid = np.concatenate([slab["resid"], first_water_resid + wat_mol])
    all_resname = list(slab["resname"]) + wat_resname
    all_atomname = list(slab["atomname"]) + wat_atomname
    all_xyz = np.vstack([slab_xyz, wat_xyz])

    _write_gro(gro_out,
               "%s | water sandwich, %.2f nm per layer" % (slab["title"].strip(), t),
               all_resid, all_resname, all_atomname, all_xyz, new_box)

    # ---------------------------------------------------------------
    # 3.8 Patch the topology
    # ---------------------------------------------------------------

    if "charmm36-jul2022" in s_forceField.lower():
        water_itp_include = "charmm36-jul2022.ff/tip3p.itp"
    elif "martini3001" in s_forceField.lower():
        pass #water in the main martini3001 itp, so no need to include a separate water itp
    elif "martini22" in s_forceField.lower():
        pass #water in the main martini22 itp, so no need to include a separate water itp
        
    
    _patch_topology(top_out, water_resname, n_water, water_itp_include)

    # ---------------------------------------------------------------
    # 3.9 Clean up and report
    # ---------------------------------------------------------------
    if not keep_workdir:
        bricksFileSystem.run_and_capture(f'rm -rf "{work_dir}"')

    if verbose:
        print("Original box : %.3f x %.3f x %.3f nm  (%d atoms)"
              % (Lx, Ly, Lz, len(slab["xyz"])))
        print("New box      : %.3f x %.3f x %.3f nm  (%d atoms)"
              % (Lx, Ly, Lz_new, len(all_xyz)))
        print("Slab now at  : z = %.3f .. %.3f nm" % (t, t + Lz))
        print("Water added  : %d molecules (%d removed for clashing with the slab)"
              % (n_water, n_mol - n_water))
        print("Wrote        : %s" % gro_out)
        print("               %s" % top_out)


    # ---------------------------------------------------------------
    # 4. bring forcefield  to the folder
    # ---------------------------------------------------------------

    if s_forceField in ["charmm36-jul2022", "martini3001", "martini22"]: #if the user chose one the forcefields that I have stored myself in USEFUL_FORCEFIELDS
    
        module_path = Path(__file__).resolve().parent # Where the cl module lives
        s_ffLocation = module_path / "USEFUL_FORCEFIELDS"
        
        if "charmm36-jul2022" in s_forceField.lower():
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/charmm36-jul2022.ff" {out_dir.rstrip("/")}/') # in the case of charmm, the actual folder has a .ff in the end
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/toppar" {out_dir.rstrip("/")}/') #and this extra file must alse come
            #ff_inclusion_text = "charmm36-jul2022.ff/forcefield.itp"
        elif "martini3001" in s_forceField.lower():
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/martini3001" {out_dir.rstrip("/")}/')
            #ff_inclusion_text = "martini3001/martini_v3.0.0.itp"
        elif "martini22" in s_forceField.lower():
            bricksFileSystem.run_and_capture(f'cp -r "{s_ffLocation}/martini22" {out_dir.rstrip("/")}/')
            #ff_inclusion_text = "martini22/martini_v2.2.itp"
        
        
    else: #if the user gave the folder of the forcefield
        s_forceField = str(Path(s_forceField).expanduser().resolve()) #resolve all ../ ~/ ../../ ./ to absolute path
        s_ffLocation = bricksFileSystem.get_file_location(s_forceField) # get forcefiled original location (relative to where the program was louched)
        s_forceField_name = bricksFileSystem.get_filename_without_extension(s_forceField) # now we update the variable so It will have just the ff name. without location nor extention
        
        bricksFileSystem.run_and_capture(f'cp -r "{s_forceField}" {out_dir.rstrip("/")}/') # in the case of charmm, the actual folder has a .ff in the end
        ff_inclusion_text = f"{s_forceField_name}/xxxxxxxx.itp"

    
 
    #copy usefull scripts to the system folder
    files_to_copy = [
        "runREALISTIC.sh",
        "runREALISTIC_expandz_cutoff.sh",
        "runREALISTIC_expandz_pme.sh",
        "runREALISTIC_nvt_cutoff.sh",
        "runREALISTIC_nvt_pme.sh",
    ]
    module_path = Path(__file__).resolve().parent # Where the cl module lives
    source_dir = module_path / "bash" # Folder containing the source files
    dest_dir = Path(out_dir) # Destination folder
    for filename in files_to_copy:
        src = source_dir / filename
        dst = dest_dir / filename
        shutil.copy(src, dst)
