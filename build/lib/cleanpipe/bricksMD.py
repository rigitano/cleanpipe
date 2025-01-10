from cleanpipe import bricksFileSystem
from cleanpipe import bricksTOP

import subprocess
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


@ensure_original_directory
def pdb2system(s_pdbfile,s_outName,s_forceField,s_boxSize,b_addterminal=True):
    """
    creates a new folder with the system name. and a gro and top files inside it with that same system name
    the top will be a socked top, all the molecules will be outside

    

    b_addterminal : True to add the standard termini and avoid dangling bonds, False to use when there are already termini in the input pdb, so no extra termini addition is required
    
    """



    #check if the filename inside s_pdbfile is valid
    bricksFileSystem.check_extention(s_pdbfile,['.pdb']) 

    #get the pdb basename. it should be the name of the protagonist molecule
    s_molName = s_pdbfile.replace(".pdb","")

    #create output folder in parael with the pdb input. we will cd into that forder and do everithing there
    bricksFileSystem.run_and_capture(f"mkdir {s_outName}")
    bricksFileSystem.run_and_capture(f"cp {s_pdbfile} {s_outName.rstrip('/')}/temp.pdb")
    bricksFileSystem.run_and_capture(f"cp -r {s_forceField}.ff {s_outName.rstrip('/')}/")#copy the forcefield to the new folder
    #original_directory = os.getcwd()#original folder is stored so I can go back to it at the very end of this function
    os.chdir(f"{s_outName}")



    ################################## create gro and top from pdb. then add the box size to the gro ###########################

    #pdb2gmx

    if b_addterminal == True:
        #standar option, that adds the correct termini in proteins
        bricksFileSystem.run_and_capture(f"gmx pdb2gmx -f temp.pdb -o {s_outName}.gro -p {s_outName}.top -i {s_molName}.posres.itp -missing -ignh -water none -ff {s_forceField}")
    elif b_addterminal == False:
        #this option will leave dangling bonds. its usefull just in case I will add termini manually
        bricksFileSystem.run_and_capture(f"printf '8\n7\n' | gmx pdb2gmx -f temp.pdb -o {s_outName}.gro -p {s_outName}.top -i {s_molName}.posres.itp -missing -ter -ignh -water none -ff {s_forceField}")
    

    #pdb2gmx is stupid, so by default it and givesa wierd name to the molecule from the pdb. most times is "Other_chain_O". lets replace it by the real molecule name, that I took from the pdb file name
    uglyMolName = bricksTOP.getMoleculeName(f"{s_outName}.top")
    bricksTOP.replaceMoleculeName(f"{s_outName}.top", uglyMolName, s_molName)

    #define box size inside the gro file. s_boxSize contains the user definition (ex: "3 3 3")
    bricksFileSystem.run_and_capture(f"gmx editconf -f {s_outName}.gro -o {s_outName}.gro -c -box {s_boxSize} -bt cubic")
    bricksFileSystem.delete(f"\\#{s_outName}.gro.1\\#")# I chose to overwrite the old gro

    #decompose the original top into a new top and a itp. the new top will contain just sytem information, the itp will describe the protagonist molecule
    bricksTOP.decompose_TOP_file_into_TOP_and_ITPs(f"{s_outName}.top")

    #delete the temporary pdb used by pdb2gmx as its no longer necessary
    bricksFileSystem.delete(f"temp.pdb")

    #after performing the system creation, go back to the original folder python was called
    #os.chdir(original_directory)

@ensure_original_directory
def solvate_and_neutralize(s_systemFolder,s_solventName,s_forceField):
    """
    there are two ways to solvate in gromacs:

        "gmx solvate -cp SoluteMolecule.gro -cs preEquilibratedBoxOfSmallSolvents.gro -o outSystem.gro -p soluteMolecule.top", 
            this will insert a molecule into a box of pre equilibrated smal solvents with just one residue. 
            after the insertion overlaping molecules will be deleted, thats why this is not a good option for large molecules of solvent, such as octane. 
            but this is a great option to solvate something into a mixture of small molecules
            I hate that the option -p soluteMolecule.top will just update the solvent molecule number in the solvent top, without including the solvent itp. the itp inclusion has to be done manualy
            if you dont put -cs, the tool will use a spc216.gro box stored in the shared/gromacs/top folder. this gro can be used to solvate any 3 other point water, such as the famous tip3p 

    
        "gmx insert-molecules -f SoluteMolecule.gro -ci solventMoleculeToInsert.gro -nmol 1000 -rot -box 5 5 5 -o outSystem.gro"
            will randomly insert solvent molecules around the solute

        s_systemFolder :system to be solvated, this mean the input is a system with only the protagonist solute that must be solvated
        s_solventName : the name of the solvent. It has to be one of gromacs standard water names, or a folder containing a preequilibrated system that is a box full of something. in both cases,from that name the function will find the necessary .gro and .itp somewere. the gro have to describe a box full of that solvent
    """
    print(f"\nCLEANPIPE MESSAGE called solvate_and_neutralize({s_systemFolder},{s_solventName},{s_forceField})\n")

    #get gro basaname in system folder
    s_groName = bricksFileSystem.get_single_gro(s_systemFolder).replace('.gro','')
    #get top basaname in system folder
    s_topName = bricksFileSystem.get_single_top(s_systemFolder).replace('.top','')


    if s_solventName in ["tip3p", "spc", "spce"]: #this is a list of 3 point water models. their respectives .gro describing a pre-equilibrated box and .itp are already in the share/gromacs/top folder
        print(f"\nCLEANPIPE MESSAGE user chose one of the standard water models ({s_solventName})\n")
        #this mean the user has chosen a water model, already part of gromacs standard solvents. gromacs can find the solvent box and the respective itp automaticaly

        #go to system folder. the current folder is savad so to go back to it just before the end of the function
        #original_directory = os.getcwd()
        os.chdir(f"{s_systemFolder}")

        bricksFileSystem.run_and_capture(f"gmx solvate -cp {s_groName}.gro -cs spc216.gro -p {s_topName}.top -o {s_groName}.gro") # spc216.gro is a pre-equilibrated box of a 3 point water model that can be used by any other 3 point model
        bricksFileSystem.delete(f"\\#{s_groName}.gro.1\\#")#I choose to overwrite the old gro
        bricksFileSystem.delete(f"\\#{s_topName}.top.1\\#")#I choose to overwrite the old top



        #include necessary text in the top file
        s_text_to_insert = "\n; Include water topology\n#include \""+s_forceField+".ff/"+s_solventName+".itp\"\n\n#ifdef POSRES_WATER \n; Position restraint for each water oxygen\n[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n1    1       1000       1000       1000\n#endif\n"

        bricksTOP.insert_text_before_directive(f"{s_topName}.top", s_text_to_insert, "[ system ]")


        # xxx add ions, to make the box neutral
        #subprocess.run(f"gmx grompp -f ~/mdparameters/add_ions.mdp -c coord_box_sol.gro -p topol.top -o coord_box_sol_ions.tpr" , shell=True, check=True)
        #subprocess.run(f"printf '13' | gmx  genion -s coord_box_sol_ions.tpr -o coord_box_sol_ions.gro -p topol.top -pname NA -nname CL -neutral" , shell=True, check=True)
        #rm tpr and created backups and mdout


    elif bricksFileSystem.check_folder(os.path.abspath(s_solventName)) == True:
        # this mean the user has chosen a folder (ex: path/to/folder)
        # that folder shoulrd contain a system that is a box filled with solvent. it should be pre-equilibrated 
        # so, the solvent name is something like box_full_of_octn, and that folder should contain a octn.itp and a 3_NPT/box_full_of_octn.gro
        # but dont worry about the gro and file names. the important is that they are present in the correct place. the name will be obtained

        #get the full path
        s_solventFolder = os.path.abspath(s_solventName)

        #obtain the names of the top and itp files in the SOLVENT BOX folder
        s_solbox_groName = bricksFileSystem.get_single_gro(f"{s_solventFolder}/3_NPT")
        l_solbox_itpNames = bricksFileSystem.get_all_itps(s_solventFolder)

        #go to system folder. the current folder is savad so to go back to it just before the end of the function
        #original_directory = os.getcwd()
        os.chdir(f"{s_systemFolder}")

        #insert the solvent in gro. and inform quantity added in top
        bricksFileSystem.run_and_capture(f"gmx solvate -cp {s_groName}.gro -cs {s_solventFolder}/3_NPT/{s_solbox_groName} -p {s_topName}.top -o {s_groName}.gro")
        bricksFileSystem.delete(f"\\#{s_groName}.gro.1\\#")#I choose to overwrite the old gro
        bricksFileSystem.delete(f"\\#{s_topName}.top.1\\#")#I choose to overwrite the old top

        #when gmx solvate inform the quantity added in the top, its possible that it chooses a weird name. lets make sure its the name of the itp file
        badmolName = bricksTOP.getMoleculeName(f"{s_topName}.top", order=-1)#get name of the last molecule in the directive [ molecules ]
        bricksTOP.replaceWordInsideDirective(f"{s_topName}.top", "[ molecules ]", badmolName, l_solbox_itpNames[0].replace(".itp", ""))# xxx this is not preparet to deal with a solvent box with several different molecules

        #edit top to insert a line including a reference of the solvent itp before the [ system ] directive
        for s_sol_itpName in l_solbox_itpNames:
            subprocess.run(rf'''awk -v line='#include "{s_sol_itpName}"' '/\[ system \]/{{print line"\n"; i=2}}i&&!--i{{next}}1' {s_topName}.top > temp.top && mv temp.top {s_topName}.top''', shell=True, check=True)

        #copy all the itp files from the original folder to the current system folder
        for s_sol_itpName in l_solbox_itpNames:
            subprocess.run(f"cp {s_solventFolder}/{s_sol_itpName} ./" , shell=True, check=True)


    else:
        print("solvation failed")





@ensure_original_directory
def make_realistic(s_systemFolder,s_groups_to_monitor_separately, s_temperature):
    """
    usage example:
    cl.make_realistic("1LZ1_in_water", s_groups_to_monitor_separately="Protein Non-Protein", s_temperature="300")
    cl.make_realistic("octn_filledbox", s_groups_to_monitor_separately="System", s_temperature="300")

    gromacs have a idiotic way to make a system realistic, you have to EM, and then simulate it 2 times, but holding the protein in place, and with several spetial mdp paramenters, so not to missfold the protein
    this is a function to do all of that without having to think about it every time
    

    will perform EM NVT and NPT in a system, so to make it realistic
    I have to give the name of a folder that contains a system
    the function will create folders 1_EM, 2_NVT and 3_NPT inside that given folder
    after this function is run, the .gro inside the 3_NPT will be the one to use in the future production runs, wereas the .top will remain the same

    s_folder : name of a folder that contains a system. this mean inside the folder there is a gro, a top, and they contain a solvated box
    s_groups_to_monitor_separately : "System" or "Protein Non-Protein"
    s_temperature: temerature in kelvin to be set during NVT, and kept during NPT, for example "300"
    
    """
    print(f"\nCLEANPIPE MESSAGE called make_realistic({s_systemFolder},{s_groups_to_monitor_separately}, {s_temperature})\n")


    #change current folder to the system folder
    os.chdir(f"{s_systemFolder}")


    # setup the correct mdp files to use, accorting to the groups, temperature, and force field
    s_module_folder = os.path.dirname(__file__)
    s_mdp_folder = os.path.join(s_module_folder,"mdp")


    #file_em_mdp="em.mdp"
    #file_nvt_mdp="nvt_begin_Vr.mdp"
    #file_npt_mdp="npt_continuation_Vr_Cr.mdp"
    #file_prod_mdp="prod_continuation_Vr_Cr.mdp"


    if s_groups_to_monitor_separately == "Protein Non-Protein":
        s_mdpNameNVT = "begin_mdNvt_Vr.mdp"
        s_mdpNameNPT = "continue_mdNpt_Vr_PaRa.mdp"
    elif s_groups_to_monitor_separately == "System":
        s_mdpNameNVT = "nvt_begin_Vr_1GROUP.mdp"
        s_mdpNameNPT = "npt_begin_Vr_Cr_1GROUP.mdp"  

    # obtain gro and top in current folder
    s_initialgroName = bricksFileSystem.get_single_gro(".")
    s_topName        = bricksFileSystem.get_single_top(".")



    #EM
    bricksFileSystem.run_and_capture(f"mkdir 1_EM")
    os.chdir(f"1_EM")
    bricksFileSystem.run_and_capture(f"gmx grompp -f {s_mdp_folder}/em.mdp -c ../{s_initialgroName} -p ../{s_topName} -o em.tpr -maxwarn 3" )
    bricksFileSystem.run_and_capture(f"gmx mdrun -deffnm em" )
    os.chdir(f"..")

    #NVT equilibration
    bricksFileSystem.run_and_capture(f"mkdir 2_NVT")
    os.chdir(f"2_NVT")
    bricksFileSystem.run_and_capture(f"gmx grompp -f {s_mdp_folder}/{s_mdpNameNVT} -c ../1_EM/em.gro -r ../1_EM/em.gro -p ../{s_topName} -o nvt.tpr -maxwarn 3")
    bricksFileSystem.run_and_capture(f"gmx mdrun -deffnm nvt")
    os.chdir(f"..")

    #NPT equilibration
    bricksFileSystem.run_and_capture(f"mkdir 3_NPT")
    os.chdir(f"3_NPT")
    bricksFileSystem.run_and_capture(f"gmx grompp -f {s_mdp_folder}/{s_mdpNameNPT} -c ../2_NVT/nvt.gro -r ../2_NVT/nvt.gro -t ../2_NVT/nvt.cpt -p ../{s_topName} -o npt.tpr -maxwarn 3")
    bricksFileSystem.run_and_capture(f"gmx mdrun -deffnm npt")
    os.chdir(f"..")



#def benchmark():

#def benchmark_rome():




@ensure_original_directory
def production(s_systemFolder,n_nanoseconds):
    """
    Runs molecular dynamics simulation for the given system folder.

    Parameters:
    s_systemFolder (str): The folder containing the system files for the simulation.
    """

    print(f"\nCLEANPIPE MESSAGE called production({s_systemFolder},{n_nanoseconds})\n")

    n_nanoseconds = int(n_nanoseconds)

    # setup mdp
    s_module_folder = os.path.dirname(__file__)
    s_mdp_folder = os.path.join(s_module_folder,"mdp")

    # go inside the system folder and create 4_PROD
    os.chdir(f"{s_systemFolder}")#change current folder to the system folder
    s_topName = bricksFileSystem.get_single_top(".") #get top name on the system folder
    bricksFileSystem.create_folder("4_PROD")#create folder
    os.chdir(f"4_PROD")#go to created folder

    #initial simulation
    bricksFileSystem.run_and_capture(f"gmx grompp -f {s_mdp_folder}/prod_continuation_Vr_Cr.mdp -c ../3_NPT/npt.gro -p ../{s_topName} -o prod.tpr -maxwarn 3")
    bricksFileSystem.run_and_capture(f"gmx mdrun -deffnm prod")

    #extend the simulation, use the tpbconv tool to extend the .tpr file
    bricksFileSystem.run_and_capture(f"gmx convert-tpr -s prod.tpr -extend {str(n_nanoseconds*1000-1)} -o prod.tpr")
    bricksFileSystem.run_and_capture(f"gmx mdrun -deffnm prod -cpi prod.cpt")  # Continue the simulation from the checkpoint file

    #xxx room for improvement in the processing configurations, for example -nt 8

    #recenter (dont go to edges) and fit (appear to just jitter standing still)
    bricksFileSystem.run_and_capture(f"printf '1\n0' | gmx trjconv -s prod.tpr -f prod.xtc -o prod_centered.xtc -center -pbc mol")
    bricksFileSystem.run_and_capture(f"printf '1\n0' | gmx trjconv -s prod.tpr -f prod_centered.xtc -o prod_fitted.xtc -fit progressive")

    
    #  fool proff alternative
    #subprocess.run(["gmx", "trjconv", "-s", "prod.tpr", "-f", "prod.xtc", "-o", "prod_centered.xtc", "-center", "-pbc", "mol"], input="1\n0\n", text=True)
    #subprocess.run(["gmx", "trjconv", "-s", "prod.tpr", "-f", "prod_centered.xtc", "-o", "prod_fitted.xtc", "-fit", "progressive"], input="1\n0\n", text=True)

    os.chdir(f"..")


#def production_rome():


#def production_fep():

#def production_fep_rome():



def hbonds(s_xtc,s_tpr):
    """
    it will run 3 hydrogen bond calculations: solvent, solvent, protein-solvent and protein-protein
    for that, a choice will have to be made in the prompt xxx but I dont know if those choices change from system to system!
    
    example
    cl.hbonds("alaEW_prod_298_00.all.fitted.xtc","alaEW_prod_298_00.tpr")
    """

    s_xtc_folder = bricksFileSystem.get_file_location(s_xtc)
    bricksFileSystem.create_folder(s_xtc_folder + '/analysis')
    s_out_folder = s_xtc_folder + '/analysis'

    

    #protein to solvent
    bricksFileSystem.run_and_capture(f"printf \"1\n13\n\" | gmx hbond -f {s_xtc} -s {s_tpr} -num {s_out_folder}/hb_ps.xvg")
    #protein to protein
    bricksFileSystem.run_and_capture(f"printf \"1\n1\n\" | gmx hbond -f {s_xtc} -s {s_tpr} -num {s_out_folder}/hb_pp.xvg")
    #solvent to solvent
    bricksFileSystem.run_and_capture(f"printf \"13\n13\n\" | gmx hbond -f {s_xtc} -s {s_tpr} -num {s_out_folder}/hb_ss.xvg")

def sasa(s_xtc,s_tpr):
    """
    
    example
    cl.sasa("alaEW_prod_298_00.all.fitted.xtc","alaEW_prod_298_00.tpr")
    """

    s_xtc_folder = bricksFileSystem.get_file_location(s_xtc)
    bricksFileSystem.create_folder(s_xtc_folder + '/analysis')
    s_out_folder = s_xtc_folder + '/analysis'

    bricksFileSystem.run_and_capture(f"printf \"1\n\" | gmx sasa -f {s_xtc} -s {s_tpr} -o {s_out_folder}/sasa.xvg")

def rama(s_xtc,s_tpr):
    """
    
    example
    cl.rama("alaEW_prod_298_00.all.fitted.xtc","alaEW_prod_298_00.tpr")
    """

    s_xtc_folder = bricksFileSystem.get_file_location(s_xtc)
    bricksFileSystem.create_folder(s_xtc_folder + '/analysis')
    s_out_folder = s_xtc_folder + '/analysis'

    bricksFileSystem.run_and_capture(f"gmx rama -f {s_xtc} -s {s_tpr} -o {s_out_folder}/rama.xvg")

    # Read the rama.xvg file and filter out lines containing @ or #, then write the first and second columns to a CSV file.
    with open(f"{s_out_folder}/rama.xvg", "r") as infile, open(f"{s_out_folder}/rama.csv", "w") as outfile:
        for line in infile:
            if not line.startswith(("@", "#")):
                columns = line.split()
                outfile.write(f"{columns[0]},{columns[1]}\n")

def dssp(s_xtc,s_gro):
    """
    
    example
    cl.dssp("alaEW_prod_298_00.all.fitted.xtc","alaEW_npt.gro")
    """

    s_xtc_folder = bricksFileSystem.get_file_location(s_xtc)
    bricksFileSystem.create_folder(s_xtc_folder + '/analysis')
    s_out_folder = s_xtc_folder + '/analysis'

    bricksFileSystem.run_and_capture(f"gmx dssp -f {s_xtc} -s {s_gro} -o {s_out_folder}/dssp.dat -hmode dssp")

    




