import subprocess
import re
from cleanpipe import filemanager
from cleanpipe import topContent



def pdb2filledBox(s_pdbfile):
    """
    #fillBox()?
    create a 5x5x5 box system filled with a lot of copies of the molecule
    """

    #check if the filename inside s_pdbfile is valid
    filemanager.check_file(s_pdbfile,['.pdb']) 
    #obtein just the file name. ex: blabla/blabla/filename.bla
    s_filename = filemanager.get_filename_without_extension(s_pdbfile) 


    #create a system with 1 molecule
    subprocess.run(f"gmx pdb2gmx -f {s_filename}.pdb -o {s_filename}.gro -p {s_filename}.top -i {s_filename}_posres.itp -water none -ff charmm36-jul2022" , shell=True)
    
    #manipulate the GRO file to create a 5x5x5 box and fill it with copyes of the molecule
    result = subprocess.run(f"gmx insert-molecules -ci {s_filename}.gro -nmol 1000 -box 5 5 5 -o {s_filename}_filledbox.gro" , shell=True, capture_output=True,text=True)
    print(result.stdout+result.stderr)

    
    #get the number of molecules realy added. this will be done by reading the standard output
    match = re.search(r'Added\s+(\d+)\s+molecules', result.stdout+result.stderr)
    added_molecules = int(match.group(1))

    #change the molecule name inside the TOP file the same as the original pdb file
    uglyMolName = topContent.getMoleculeName({s_filename}+".top")
    molName = s_filename
    topContent.replaceMoleculeName({s_filename}+".top", uglyMolName, molName)

    #update the TOP file with the new total the molecule
    topContent.update_molecule_quantity({s_filename}+".top", molName, added_molecules)


