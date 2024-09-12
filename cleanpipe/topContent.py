import os
import subprocess
import re

def getMoleculeName(top_file_path, order=1):
    """
    obtains the name of a molecule inside a top file. 
    You have to give the top file
    You can give the molecule order, if there are several molecules in the top file. To get the last, give -1. 
    """
    molecules = []
    
    # Open and read the .top file
    with open(top_file_path, 'r') as file:
        lines = file.readlines()
        
        # Flag to indicate if we are in the [ molecules ] section
        in_molecules_section = False
        
        for line in lines:
            line = line.strip()
            
            # Check if the line is the start of the [ molecules ] section
            if line.startswith('[ molecules ]'):
                in_molecules_section = True
                continue
            
            # If we are in the [ molecules ] section and encounter a non-comment, non-empty line
            if in_molecules_section and line and not line.startswith(';'):
                # Extract the first word in the line, which is the molecule name
                molecule_name = line.split()[0]
                molecules.append(molecule_name)
    
    # Handle the case for choosing the last molecule if order is -1
    if order == -1:
        return molecules[-1] if molecules else None
    
    # Handle the case for other order values (1-based index)
    index = order - 1
    if 0 <= index < len(molecules):
        return molecules[index]
    else:
        return None  # If the requested order is out of bounds
    
def getSystemName(top_file_path):
    """
    obtains the name of the system inside a top file. 
    You have to give the top file

    """
    
    # Open and read the .top file
    with open(top_file_path, 'r') as file:
        lines = file.readlines()
        
        # Flag to indicate if we are in the [ system ] section
        in_system_section = False
        
        for line in lines:
            line = line.strip()
            
            # Check if the line is the start of the [ molecules ] section
            if line.startswith('[ system ]'):
                in_system_section = True
                continue
            
            # If we are in the [ system ] section and encounter a non-comment, non-empty line
            if in_system_section and line and not line.startswith(';'):
                # Extract the first line that is not a comment, this should be the system name
                return line
    



def replaceWordInsideDirective(top_filename, target_directive, old_word, new_word):
    """
    this replaces a word inside a specific directive. 
    you gave to informe the top file, the directive to look into, the wor to be replaced, and the new word
    """

    # Reading the file contents
    with open(top_filename, 'r') as file:
        lines = file.readlines()

    inside_target_directive = False
    updated_lines = []

    for line in lines:
        stripped_line = line.strip()

        # Look for target_directive (ex: [ moleculetype ]) and its contents
        if stripped_line.startswith(target_directive):
            inside_target_directive = True
            updated_lines.append(line)
            continue

        # Exiting target_directive (ex: [ moleculetype ]) when a blank line or another directive is found
        if inside_target_directive and stripped_line.startswith("[") and target_directive not in stripped_line:
            inside_target_directive = False

        # Replace the molecule name in the target_directive (ex: [ moleculetype ]) 
        if inside_target_directive and old_word in stripped_line:
            updated_lines.append(line.replace(old_word, new_word))
        else:
            updated_lines.append(line)

    # Writing the modified content back to the file
    with open(top_filename, 'w') as file:
        file.writelines(updated_lines)

def replaceMoleculeName(top_filename, old_molecule_name, new_molecule_name):
    """
    to replace a molecule name, you have to do replace the name that appears inside two specific directives: [ moleculetype ] and [ molecules ]
    """

    replaceWordInsideDirective(top_filename, "[ moleculetype ]", old_molecule_name, new_molecule_name)
    replaceWordInsideDirective(top_filename, "[ molecules ]", old_molecule_name, new_molecule_name)


def setSystemName(top_filename, new_system_name):
    """
    to replace the system name, this function find the current system name and then replace the name that appears inside two specific directives: [ moleculetype ] and [ molecules ]
    """

    old_system_name = getSystemName(top_filename)

    replaceWordInsideDirective(top_filename, "[ system ]", old_system_name, new_system_name)


def update_molecule_quantity(top_file, molecule_name, new_quantity):
    """
    Update the quantity of a specific molecule in the [ molecules ] directive in a GROMACS top file.

    :param top_file: Path to the GROMACS .top file
    :param molecule_name: Name of the molecule whose quantity needs to be changed
    :param new_quantity: New quantity to replace the old one
    """
    # Read the content of the top file
    with open(top_file, 'r') as file:
        lines = file.readlines()

    # Flags to indicate that we are inside the [ molecules ] directive
    in_molecules_section = False
    new_lines = []

    # Process each line
    for line in lines:
        # Check if we are entering the [ molecules ] section
        if '[ molecules ]' in line:
            in_molecules_section = True
            new_lines.append(line)
            continue

        # Exit the [ molecules ] section if we hit a blank line or a section header
        if in_molecules_section and (line.strip() == '' or line.strip().startswith('[')):
            in_molecules_section = False

        # If in the molecules section, modify the specific molecule quantity
        if in_molecules_section:
            # Split the line into the molecule name and its quantity
            split_line = line.split()
            if len(split_line) == 2 and split_line[0] == molecule_name:
                # Replace the quantity with the new one
                new_line = f"{molecule_name}    {new_quantity}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write the modified lines back to the file
    with open(top_file, 'w') as file:
        file.writelines(new_lines)


def decompose_TOP_file_into_SOCKETTOP_and_ITPs(top_file_path):

    # Read the content of the original top file
    with open(top_file_path, 'r') as f:
        lines = f.readlines()

    # Storage for system information and molecule-specific sections
    system_info = []
    molecule_sections = {}
    molecule_names = {}
    current_molecule_name = None
    current_molecule_id = -1
    inside_molecule = False
    inside_moleculetype_directive = False

    # Go through file lines
    for line in lines:

        # Check if we are in the moleculetype title
        if line.startswith("[ moleculetype ]"):
            # Start of a new molecule type
            inside_molecule = True
            inside_moleculetype_directive = True
            current_molecule_name = None  # Reset current molecule
            current_molecule_id += 1
            molecule_sections[current_molecule_id] = [] #initialise place were molecule infos are going to be stored

        #check if the molecules are all passed and we are in the system description at the end of the file
        elif line.startswith("[ system ]") or line.startswith("[ molecules ]"):
            inside_molecule = False
            inside_moleculetype_directive = False
            # Detect global system-related sections after molecule definitions

        #check if we are in any other title
        elif inside_molecule and (line.startswith("[") and not line.startswith("[ moleculetype ]")):
            # You are in the title of a directive other than [ moleculetype ]
            inside_moleculetype_directive = False

        # check if we are in the line with the molecule name
        elif inside_moleculetype_directive and line and not line.startswith(';') and bool(line.strip()): 
            current_molecule_name = line.strip().split()[0] # Extract the first word in the line, which is the molecule name
            molecule_names[current_molecule_id] = current_molecule_name

        
        if inside_molecule:
            # Append lines related to the current molecule
            molecule_sections[current_molecule_id].append(line) #store line with molecule info
        elif not inside_molecule:
            #this will be true before finding the first [ moleculetype ] directive
            system_info.append(line)
    
    # Create new files based on the parsed data
    top_dir = os.path.dirname(top_file_path)
    base_name = os.path.splitext(os.path.basename(top_file_path))[0]
    
    # Create the new system top file without molecule definitions
    system_top_file = os.path.join(top_dir, f"{base_name}.socket.top")
    with open(system_top_file, 'w') as f:
        for line in system_info:
            f.write(line)

    # Add #include for each molecule itp file
    for cont in range(len(molecule_names) - 1, -1, -1):
        molecule_name = molecule_names[cont]
        subprocess.run(rf'''awk -v line='#include "{molecule_name}.itp"' '/\[ system \]/{{print line"\n"; i=2}}i&&!--i{{next}}1' {system_top_file} > temp.top && mv temp.top {system_top_file}''', shell=True, check=True)

    print(f"top file written: \n                      {system_top_file}")
    
    # Create separate itp files for each molecule
    print(f"{len(molecule_names)} itp file written:")
    for molecule_id, section_lines in molecule_sections.items():
        itp_file = os.path.join(top_dir, f"{molecule_names[molecule_id]}.itp")
        with open(itp_file, 'w') as f:
            # Write molecule-specific content to the itp file
            f.writelines(section_lines)
        print(f"                      {itp_file}")




def remove_posres_inclusion(s_topfile):
    """
    unfortunately pdb2gmx always create a posres.itp file. this is sometimes useless, for example
    when creating a box full of the same molecule to use it as a custom solvent
    this function removes the inclusion of a posres.itp file 
    """

    print(s_topfile)

    # Create a more flexible regex pattern to match the inclusion block
    pattern_to_remove = r';\s*Include\s*Position\s*restraint\s*file\s*\n#ifdef\s*POSRES\s*\n#include\s*"posres\.itp"\s*\n#endif\s*'

    # Open the file and read its contents
    with open(s_topfile, 'r') as file:
        file_contents = file.read()
        print(file_contents)

    # Remove the matched block using the regex pattern
    updated_contents = re.sub(pattern_to_remove, '', file_contents, flags=re.MULTILINE)
    print(updated_contents)

    # Write the updated contents back to the file
    with open(s_topfile, 'w') as file:
        file.write(updated_contents)

    
     