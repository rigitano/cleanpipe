import subprocess
import re
from cleanpipe import filemanager
from cleanpipe import topContent



def pdb2filledBox(s_pdbfile):
    """
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


    #create a system with 1 molecule.
    subprocess.run(f"gmx pdb2gmx -f {s_filename}.pdb -o {s_filename}.gro -p {s_filename}.top -i {s_filename}.posres.itp -water none -ff charmm36-jul2022" , shell=True)
    
    #manipulate the GRO file to create a 5x5x5 box and fill it with copyes of the molecule
    result = subprocess.run(f"gmx insert-molecules -ci {s_filename}.gro -nmol 1000 -box 5 5 5 -o {s_filename}_filledbox.gro" , shell=True, capture_output=True,text=True)
    print(result.stdout+result.stderr)
    subprocess.run(f"rm {s_filename}.gro" , shell=True)# now that we have the filled box gro, the 1 molecule gro can be deleted

    #get the number of molecules really added. this will be done by reading the standard output
    match = re.search(r'Added\s+(\d+)\s+molecules', result.stdout+result.stderr)
    added_molecules = int(match.group(1))

    #rename the top and itp files. but notice the top will have to be edited so to reflect the new molecule total. the name of the system and molecules will also be edited
    subprocess.run(f"mv {s_filename}.top {s_filename}_filledbox.top" , shell=True)
    subprocess.run(f"mv {s_filename}.posres.itp {s_filename}_filledbox.posres.itp" , shell=True)


    #change the ugly molecule name currently inside the TOP file. it will be changed to be the same as the original pdb file
    uglyMolName = topContent.getMoleculeName(f"{s_filename}_filledbox.top")
    molName = s_filename
    topContent.replaceMoleculeName(f"{s_filename}_filledbox.top", uglyMolName, molName)
    topContent.replaceMoleculeName(f"{s_filename}_filledbox.top", uglyMolName, molName)

    #update the TOP file with the new total the molecule
    topContent.update_molecule_quantity(f"{s_filename}_filledbox.top", molName, added_molecules)

    #give a name for the system
    topContent.setSystemName(f"{s_filename}_filledbox.top", f"box filled with {s_filename}" )


