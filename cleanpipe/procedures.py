import subprocess
import re
from cleanpipe import filemanager
from cleanpipe import topContent



def pdb2filledBox(s_pdbfile, s_forceField):
    """
    usage example:
    cl.pdb2filledBox("octn.pdb","charmm36-jul2022")

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


    subprocess.run(f"mkdir {s_filename}_filledbox", shell=True)
    s_outPathAndName = f"{s_filename}_filledbox/{s_filename}"


    #create a system with 1 molecule.
    subprocess.run(f"gmx pdb2gmx -f {s_filename}.pdb -o {s_outPathAndName}.gro -p {s_outPathAndName}.top -i posres.itp -water none -ff {s_forceField}" , shell=True)
    
    #pdb2gmx generates a useless posres.itp with useless posres for 1 molecule. so I delete the posres.itp and the inclusion in the top
    subprocess.run(f"rm posres.itp" , shell=True) 
    topContent.remove_posres_inclusion(f"{s_outPathAndName}.top")

    #manipulate the GRO file to create a 5x5x5 box and fill it with copyes of the molecule
    result = subprocess.run(f"gmx insert-molecules -ci {s_outPathAndName}.gro -nmol 1000 -box 5 5 5 -o {s_outPathAndName}_filledbox.gro" , shell=True, capture_output=True,text=True)
    print(result.stdout+result.stderr)
    subprocess.run(f"rm {s_outPathAndName}.gro" , shell=True)# now that we have the filled box gro, the 1 molecule gro can be deleted
    print(f"gro file written: \n                      {s_outPathAndName}_filledbox.gro")

    #get the number of molecules really added. this will be done by reading the standard output
    match = re.search(r'Added\s+(\d+)\s+molecules', result.stdout+result.stderr)
    added_molecules = int(match.group(1))

    #rename the top and posres.itp files. but notice the top will have to be edited so to reflect the new molecule total. the name of the system and molecules will also be edited
    subprocess.run(f"mv {s_outPathAndName}.top {s_outPathAndName}_filledbox.top" , shell=True)

    #change the ugly molecule name currently inside the TOP file. it will be changed to be the same as the original pdb file
    uglyMolName = topContent.getMoleculeName(f"{s_outPathAndName}_filledbox.top")
    molName = s_filename
    topContent.replaceMoleculeName(f"{s_outPathAndName}_filledbox.top", uglyMolName, molName)
    topContent.replaceMoleculeName(f"{s_outPathAndName}_filledbox.top", uglyMolName, molName)

    #update the TOP file with the new total the molecule
    topContent.update_molecule_quantity(f"{s_outPathAndName}_filledbox.top", molName, added_molecules)

    #split the TOP file, into a ITP that describes the molecule and a simple TOP that contains only name of the system and the totals.
    topContent.decompose_TOP_file_into_SOCKETTOP_and_ITPs(f"{s_outPathAndName}_filledbox.top")
    subprocess.run(f"rm {s_outPathAndName}_filledbox.top" , shell=True)

    #give a name for the system
    topContent.setSystemName(f"{s_outPathAndName}_filledbox.socket.top", f"box filled with {s_filename}" )


def pdb2molecule_in_solvent(s_pdbfile, s_outName, s_solvent, s_forceField, s_boxSize):
    """
    usage example:
    cl.pdb2molecule_in_solvent("octn.pdb", "octane_in_solvent", "tip3p", "charmm36-jul2022", "3 3 3")
    

    s_pdbfile       : string with the pdb name. for example "insulin.pdb", this will be the main molecule in the system.
    s_outFileName   : string with the name of the system, for example "alaHW". a folder with that name will be created, and inside it, all the files, for example: alaHW.gro and alaHW.top
    s_solvent       : choose a water model, for example as "tip3p", or a folder, for example "octn_filledbox". The folder have to contain a system with a solvent box, in other words, it has to contain a octn_filledbox.gro and a octn.itp
    s_forceField    : one of the gromacs recognized force fields, for example "charmm36-jul2022"
    s_boxSize       : string with x y z sizes, for example "3 3 3"
    """

     #check if the filename inside s_pdbfile is valid
    filemanager.check_file(s_pdbfile,['.pdb']) 


    subprocess.run(f"mkdir {s_outName}", shell=True)
    s_outPathAndName = f"{s_outName}/{s_outName}"

    ################################## create gro and top from pdb. then add the box size to the gro ###########################

    #pdb2gmx
    subprocess.run(f"printf '8\n7\n' | gmx pdb2gmx -f {s_pdbfile} -o {s_outPathAndName}.gro -p {s_outPathAndName}.top -i {s_outName}.posres.itp -missing -ter -ignh -water tip3p -ff {s_forceField}", shell=True)
    #pdb2gmx generates a posres.itp and put a #include statement it in the top. the generation must be done outside the out folder so not to mess up the reference in the #include statamente
    subprocess.run(f"mv {s_outName}.posres.itp {s_outPathAndName}.posres.itp" , shell=True) 
    
    #define box size. s_boxSize contains the user definition (ex: "3 3 3")
    subprocess.run(f"gmx editconf -f {s_outPathAndName}.gro -o {s_outPathAndName}.gro -c -box {s_boxSize} -bt cubic", shell=True)
    subprocess.run(f"rm {s_outName}/\\#{s_outName}.gro.1\\#" , shell=True)# I chose to overwrite the old gro


    ################################## add solvent to the system. I have 2 options here: tip3p or filled box ##############################################

    if s_solvent == "tip3p":
        #this mean the user has chosen the tip3p water model, already part of the forcefield

        subprocess.run(f"gmx solvate -cp {s_outPathAndName}.gro -p {s_outPathAndName}.top -o {s_outPathAndName}.gro", shell=True)
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.gro.1\\#" , shell=True)#I chose to overwrite the old gro
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.top.1\\#" , shell=True)#I chose to overwrite the old top

        # xxx add ions, to make the box neutral
        #subprocess.run(f"gmx grompp -f ~/mdparameters/add_ions.mdp -c coord_box_sol.gro -p topol.top -o coord_box_sol_ions.tpr" , shell=True)
        #subprocess.run(f"printf '13' | gmx  genion -s coord_box_sol_ions.tpr -o coord_box_sol_ions.gro -p topol.top -pname NA -nname CL -neutral" , shell=True)
        #rm tpr and created backups

    elif filemanager.check_folder(s_solvent) == True:
        # this mean the user has chosen a folder (ex: path/to/folder/octn)
        # lets hope that folder contains a system that is a box filled with solvent. for example octn.gro and octn.itp


        #get the names of the top and itp files in the solvent box folder
        s_solbox_groName = filemanager.get_single_gro(s_solvent)
        l_itpNames = filemanager.get_all_itps(s_solvent)

        #edit gro to insert the solvent
        subprocess.run(f"gmx solvate -cp {s_outPathAndName}.gro -cs {s_solvent}/{s_solbox_groName} -p {s_outPathAndName}.top -o {s_outPathAndName}.gro", shell=True)
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.gro.1\\#" , shell=True)#I chose to overwrite the old gro
        subprocess.run(f"rm {s_outName}/\\#{s_outName}.top.1\\#" , shell=True)#I chose to overwrite the old top

        #edit top to insert a line including a reference of the solvent itp before the [ system ] directive
        for s_sol_itpName in l_itpNames:
            subprocess.run(rf'''awk -v line='#include "{s_sol_itpName}"' '/\[ system \]/{{print line"\n"; i=2}}i&&!--i{{next}}1' {s_outPathAndName}.top > temp.top && mv temp.top {s_outPathAndName}.top''', shell=True, check=True)

        #copy all the itp files from the original folder to the current system folder
        for s_sol_itpName in l_itpNames:
            subprocess.run(f"cp {s_solvent}/{s_sol_itpName} {s_outName}/" , shell=True)

    else:
        print("solvation failed")





    ################################### set the the name of the system in the top file #######################################
    topContent.setSystemName(f"{s_outPathAndName}.top", f"{s_outName} (molecule from {s_pdbfile}, inserted in solution made using {s_solvent})" )


def download_and_clean_pdb(s_molecule_name):
    """
    usage example:
    cl.download_and_clean_pdb("1aki")

    
    """

    #get pdb from portal
    subprocess.run(f"wget https://files.rcsb.org/download/{s_molecule_name}.pdb" , shell=True)

    #remove water
    subprocess.run(f"grep -v 'HOH' {s_molecule_name}.pdb > {s_molecule_name}_temp.pdb" , shell=True)
    subprocess.run(f"rm {s_molecule_name}.pdb" , shell=True)
    subprocess.run(f"mv {s_molecule_name}_temp.pdb {s_molecule_name}.pdb" , shell=True)


def make_realistic(s_folder):
    """
    usage example:
    cl.make_realistic("insulin_in_water")

    will perform EM NVT and NPT in a system, so to make it realistic
    I have to give the name of a folder that contains a system
    the function will create folders 1_EM, 2_NVT and 3_NPT inside that given folder
    after this function is run, the .gro inside the 3_NPT will be the one to use in the future production runs, wereas the .top will remain the same

    s_folder : name of a folder that contains a system. this mean inside the folder there is a gro, a top, and they contain a solvated box

    
    """

    s_initialgroName = filemanager.get_single_gro(s_folder)
    s_topName        = filemanager.get_single_top(s_folder)

    #subprocess.run(f"xxxxxxx" , shell=True)

    #EM
    subprocess.run(f"mkdir {s_folder}/1_EM" , shell=True)
    subprocess.run(f"gmx grompp -f ~/mdp/em.mdp -c {s_folder}/{s_initialgroName} -p {s_folder}/{s_topName} -o {s_folder}/1_EM/em.tpr -maxwarn 3" , shell=True)
    subprocess.run(f"gmx mdrun -deffnm {s_folder}/1_EM/em" , shell=True)

    #NVT equilibration
    subprocess.run(f"mkdir {s_folder}/2_NVT" , shell=True)
    subprocess.run(f"gmx grompp -f ~/mdp/begin_mdNvt_Vr.mdp -c {s_folder}/1_EM/em.gro -r  {s_folder}/1_EM/em.gro -p {s_folder}/{s_topName} -o {s_folder}/2_NVT/nvt.tpr" , shell=True)
    subprocess.run(f"gmx mdrun -deffnm {s_folder}/2_NVT/nvt" , shell=True)

    #NPT equilibration
    subprocess.run(f"mkdir {s_folder}/3_NPT" , shell=True)
    subprocess.run(f"grompp -f ~/mdp/continue_mdNpt_Vr_PaRa.mdp -c {s_folder}/2_NVT/nvt.gro -r {s_folder}/2_NVT/nvt.gro -t {s_folder}/2_NVT/nvt.cpt -p {s_folder}/{s_topName} -o {s_folder}/3_NPT/npt.tpr" , shell=True)
    subprocess.run(f"mdrun -deffnm {s_folder}/3_NPT/npt" , shell=True)



