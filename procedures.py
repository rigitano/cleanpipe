import subprocess
import os
import filemanager



def pdb2filledBox(s_pdbfile):
    """
    createas a system that is a 5x5x5 box filled with a lot of copies of the input pdb molecule
    """

    
    filemanager.check_file(s_pdbfile,['.pdb']) #check if the filename inside s_pdbfile is valid
    s_filename = filemanager.get_filename_without_extension(s_pdbfile) #obtein just the file name. ex: blabla/blabla/filename.bla

    subprocess.run(f"gmx pdb2gmx -f {s_filename}.pdb -o {s_filename}.gro -p {s_filename}.top -i {s_filename}_posres.itp -water none -ff charmm36-jul2022" , shell=True)
    subprocess.run(f"gmx insert-molecules -ci {s_filename}.gro -nmol 1000 -box 5 5 5 -o {s_filename}_filledbox.gro" , shell=True)

