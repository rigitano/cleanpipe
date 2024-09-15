from cleanpipe import filemanager
from cleanpipe import topContent

import subprocess
import os
import functools

def ensure_original_directory(func):
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
def better_pdb2gmx(s_pdbfile,s_outName,s_forceField,s_boxSize):
    """
    creates a new folder with the system name. and a gro and top files inside it with that same system name
    the top will be a socked top, all the molecules will be outside


    
    """



    #check if the filename inside s_pdbfile is valid
    filemanager.check_file(s_pdbfile,['.pdb']) 

    #get the pdb basename. it should be the name of the protagonist molecule
    s_molName = s_pdbfile.replace(".pdb","")

    #create output folder in parael with the pdb input. we will cd into that forder and do everithing there
    subprocess.run(f"mkdir {s_outName}", shell=True, check=True)
    subprocess.run(f"cp {s_pdbfile} {s_outName}/temp.pdb", shell=True, check=True)
    #original_directory = os.getcwd()#original folder is stored so I can go back to it at the very end of this function
    os.chdir(f"{s_outName}")



    ################################## create gro and top from pdb. then add the box size to the gro ###########################

    #pdb2gmx
    subprocess.run(f"printf '8\n7\n' | gmx pdb2gmx -f temp.pdb -o {s_outName}.gro -p {s_outName}.top -i {s_molName}.posres.itp -missing -ter -ignh -water none -ff {s_forceField}", shell=True, check=True)
    
    #pdb2gmx is stupid, so by default it and givesa wierd name to the molecule from the pdb. most times is "Other_chain_O". lets replace it by the real molecule name, that I took from the pdb file name
    uglyMolName = topContent.getMoleculeName(f"{s_outName}.top")
    topContent.replaceMoleculeName(f"{s_outName}.top", uglyMolName, s_molName)

    #define box size inside the gro file. s_boxSize contains the user definition (ex: "3 3 3")
    subprocess.run(f"gmx editconf -f {s_outName}.gro -o {s_outName}.gro -c -box {s_boxSize} -bt cubic", shell=True, check=True)
    subprocess.run(f"rm \\#{s_outName}.gro.1\\#" , shell=True, check=True)# I chose to overwrite the old gro

    #decompose the original top into a new top and a itp. the new top will contain just sytem information, the itp will describe the protagonist molecule
    topContent.decompose_TOP_file_into_SOCKETTOP_and_ITPs(f"{s_outName}.top")

    #delete the temporary pdb used by pdb2gmx as its no longer necessary
    subprocess.run(f"rm temp.pdb", shell=True, check=True)

    #after performing the system creation, go back to the original folder python was called
    #os.chdir(original_directory)

@ensure_original_directory
def better_solvate(s_systemFolder,s_solvent):
    """
    
    """

    #get gro basaname in system folder
    s_groName = filemanager.get_single_gro(s_systemFolder).replace('.gro','')
    #get top basaname in system folder
    s_topName = filemanager.get_single_top(s_systemFolder).replace('.top','')


    if s_solvent == "tip3p":
        #this mean the user has chosen the tip3p water model, already part of the forcefield

        #go to system folder. the current folder is savad so to go back to it just before the end of the function
        #original_directory = os.getcwd()
        os.chdir(f"{s_systemFolder}")

        subprocess.run(f"gmx solvate -cp {s_groName}.gro -p {s_topName}.top -o {s_groName}.gro", shell=True, check=True)
        subprocess.run(f"rm \\#{s_groName}.gro.1\\#" , shell=True, check=True)#I chose to overwrite the old gro
        subprocess.run(f"rm \\#{s_topName}.top.1\\#" , shell=True, check=True)#I chose to overwrite the old top


        #include necessary text in the top file

        s_text_to_insert ="""\n; Include water topology\n#include "./charmm36-jul2022.ff/tip3p.itp"\n\n#ifdef POSRES_WATER\n; Position restraint for each water oxygen\n[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n1    1       1000       1000       1000\n#endif\n\n"""

        topContent.insert_text_before_directive(f"{s_topName}.top", s_text_to_insert, "[ system ]")


        # xxx add ions, to make the box neutral
        #subprocess.run(f"gmx grompp -f ~/mdparameters/add_ions.mdp -c coord_box_sol.gro -p topol.top -o coord_box_sol_ions.tpr" , shell=True, check=True)
        #subprocess.run(f"printf '13' | gmx  genion -s coord_box_sol_ions.tpr -o coord_box_sol_ions.gro -p topol.top -pname NA -nname CL -neutral" , shell=True, check=True)
        #rm tpr and created backups and mdout

        #after performing the solvation, go back to the original folder python was called
        #os.chdir(original_directory)

    elif filemanager.check_folder(full_path = os.path.abspath(s_solvent)) == True:
        # this mean the user has chosen a folder (ex: path/to/folder/octn)
        # that folder shoulrd contain a system that is a box filled with solvent. for example octn.gro and octn.itp

        #get the full path
        s_solventFolder = os.path.abspath(s_solvent)

        #get the names of the top and itp files in the SOLVENT BOX folder
        s_solbox_groName = filemanager.get_single_gro(s_solventFolder)
        l_solbox_itpNames = filemanager.get_all_itps(s_solventFolder)

        #go to system folder. the current folder is savad so to go back to it just before the end of the function
        #original_directory = os.getcwd()
        os.chdir(f"{s_systemFolder}")

        #insert the solvent in gro. and inform quantity added in top
        subprocess.run(f"gmx solvate -cp {s_groName}.gro -cs {s_solventFolder}/{s_solbox_groName} -p {s_topName}.top -o {s_groName}.gro", shell=True, check=True)
        subprocess.run(f"rm \\#{s_groName}.gro.1\\#" , shell=True, check=True)#I chose to overwrite the old gro
        subprocess.run(f"rm \\#{s_topName}.top.1\\#" , shell=True, check=True)#I chose to overwrite the old top

        #when gmx solvate inform the quantity added in the top, its possible that it chooses a weird name. lets make sure its the name of the itp file
        badmolName = topContent.getMoleculeName(f"{s_topName}.top", order=-1)#get name of the last molecule in the directive [ molecules ]
        topContent.replaceWordInsideDirective(f"{s_topName}.top", "[ molecules ]", badmolName, l_solbox_itpNames[0].replace(".itp", ""))# xxx this is not preparet to deal with a solvent box with several different molecules

        #edit top to insert a line including a reference of the solvent itp before the [ system ] directive
        for s_sol_itpName in l_solbox_itpNames:
            subprocess.run(rf'''awk -v line='#include "{s_sol_itpName}"' '/\[ system \]/{{print line"\n"; i=2}}i&&!--i{{next}}1' {s_topName}.top > temp.top && mv temp.top {s_topName}.top''', shell=True, check=True)

        #copy all the itp files from the original folder to the current system folder
        for s_sol_itpName in l_solbox_itpNames:
            subprocess.run(f"cp {s_solventFolder}/{s_sol_itpName} ./" , shell=True, check=True)

        #after performing the solvation, go back to the original folder python was called
        #os.chdir(original_directory)

    else:
        print("solvation failed")






def make_realistic(s_systemFolder):
    """
    usage example:
    cl.make_realistic("insulin_in_water")

    will perform EM NVT and NPT in a system, so to make it realistic
    I have to give the name of a folder that contains a system
    the function will create folders 1_EM, 2_NVT and 3_NPT inside that given folder
    after this function is run, the .gro inside the 3_NPT will be the one to use in the future production runs, wereas the .top will remain the same

    s_folder : name of a folder that contains a system. this mean inside the folder there is a gro, a top, and they contain a solvated box

    
    """

    # xxx o mdp menciona proteina e o arquivo q ser otimizado menciona outro nome para a molecula protagonista
    # xxx a pasta charm36 nao esta sendo encontrada na pasta local. preciso fazer que ela fique na pasta top sendo vista sempre

    s_initialgroName = filemanager.get_single_gro(s_systemFolder)
    s_topName        = filemanager.get_single_top(s_systemFolder)

    #subprocess.run(f"xxxxxxx" , shell=True, check=True)

    #EM
    subprocess.run(f"mkdir {s_systemFolder}/1_EM" , shell=True, check=True)
    subprocess.run(f"gmx grompp -f ~/mdp/em.mdp -c {s_systemFolder}/{s_initialgroName} -p {s_systemFolder}/{s_topName} -o {s_systemFolder}/1_EM/em.tpr -maxwarn 3" , shell=True, check=True)
    subprocess.run(f"gmx mdrun -deffnm {s_systemFolder}/1_EM/em" , shell=True, check=True)

    #NVT equilibration
    subprocess.run(f"mkdir {s_systemFolder}/2_NVT" , shell=True, check=True)
    subprocess.run(f"gmx grompp -f ~/mdp/begin_mdNvt_Vr.mdp -c {s_systemFolder}/1_EM/em.gro -r  {s_systemFolder}/1_EM/em.gro -p {s_systemFolder}/{s_topName} -o {s_systemFolder}/2_NVT/nvt.tpr" , shell=True, check=True)
    subprocess.run(f"gmx mdrun -deffnm {s_systemFolder}/2_NVT/nvt" , shell=True, check=True)

    #NPT equilibration
    subprocess.run(f"mkdir {s_systemFolder}/3_NPT" , shell=True, check=True)
    subprocess.run(f"grompp -f ~/mdp/continue_mdNpt_Vr_PaRa.mdp -c {s_systemFolder}/2_NVT/nvt.gro -r {s_systemFolder}/2_NVT/nvt.gro -t {s_systemFolder}/2_NVT/nvt.cpt -p {s_systemFolder}/{s_topName} -o {s_systemFolder}/3_NPT/npt.tpr" , shell=True, check=True)
    subprocess.run(f"mdrun -deffnm {s_systemFolder}/3_NPT/npt" , shell=True, check=True)