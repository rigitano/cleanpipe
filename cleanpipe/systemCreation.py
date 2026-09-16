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
        '#include "martini3001/martini_v3.0.0_sterols_v1.0.itp"\n'
        '\n'
    )

    top_file.write_text(
        martini_includes + top_file.read_text()
    )