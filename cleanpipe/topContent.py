

def replaceMoleculeName(top_filename, old_molecule_name, new_molecule_name):
    """
    this replaces the name of a molecule present in the top file.
    you gave to informe the top file, the name of the molecule to be replaced, and the new name
    """

    # Reading the file contents
    with open(top_filename, 'r') as file:
        lines = file.readlines()

    inside_moleculetype = False
    inside_molecules = False
    updated_lines = []

    for line in lines:
        stripped_line = line.strip()

        # Look for [ moleculetype ] directive and its contents
        if stripped_line.startswith("[ moleculetype ]"):
            inside_moleculetype = True
            updated_lines.append(line)
            continue

        # Exiting [ moleculetype ] when a blank line or another directive is found
        if inside_moleculetype and stripped_line.startswith("[") and "moleculetype" not in stripped_line:
            inside_moleculetype = False

        # Replace the molecule name in the [ moleculetype ] directive
        if inside_moleculetype and old_molecule_name in stripped_line:
            updated_lines.append(line.replace(old_molecule_name, new_molecule_name))
        else:
            updated_lines.append(line)

        # Look for [ molecules ] directive and its contents
        if stripped_line.startswith("[ molecules ]"):
            inside_molecules = True
            updated_lines.append(line)
            continue

        # Exiting [ molecules ] when a blank line or another directive is found
        if inside_molecules and stripped_line.startswith("[") and "molecules" not in stripped_line:
            inside_molecules = False

        # Replace the molecule name in the [ molecules ] directive
        if inside_molecules and old_molecule_name in stripped_line:
            updated_lines.append(line.replace(old_molecule_name, new_molecule_name))
        else:
            if not inside_moleculetype:
                updated_lines.append(line)

    # Writing the modified content back to the file
    with open(top_filename, 'w') as file:
        file.writelines(updated_lines)

    print(f"Replaced all occurrences of '{old_molecule_name}' with '{new_molecule_name}' in the file '{top_filename}'.")


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