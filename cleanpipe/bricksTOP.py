import os
import subprocess
import re
import tempfile
import sys
from pathlib import Path
from cleanpipe import bricksFileSystem
from cleanpipe import bricksStorage
from cleanpipe import lltools
from cleanpipe import bricksGRO



def get_forcefield_name(top_file_path):
    """
    Extracts the forcefield folder name from a GROMACS topology file.
    Example outputs: 'oplsaa.ff', 'charmm36-jul2022.ff', 'amber99sb-ildn.ff'
    """
    ff_pattern = re.compile(r'#include\s+"(?:\.\/)?([^/]+\.ff)/forcefield\.itp"')

    with open(top_file_path, 'r') as f:
        for line in f:
            match = ff_pattern.search(line)
            if match:
                return match.group(1)  # the forcefield folder name

    return None

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
                new_line = f"{molecule_name}    {str(new_quantity)}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write the modified lines back to the file
    with open(top_file, 'w') as file:
        file.writelines(new_lines)






def decompose_TOP_file_into_TOP_and_ITPs(top_file_path):
    """
    xxx duplicado?
    """

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

    #delete the original top file
    bricksFileSystem.delete(f"{top_file_path}")
    
    # Create the new system top file without molecule definitions
    system_top_file = os.path.join(top_dir, f"{base_name}.top")
    with open(system_top_file, 'w') as f:
        for line in system_info:
            f.write(line)

    # Add #include for each molecule itp file
    for cont in range(len(molecule_names) - 1, -1, -1):
        molecule_name = molecule_names[cont]
        subprocess.run(rf'''awk -v line='#include "{molecule_name}.itp"' '/\[ system \]/{{print line"\n"; i=2}}i&&!--i{{next}}1' {system_top_file} > temp.top && mv temp.top {system_top_file}''', shell=True, check=True)

    print(f"\nCLEANPIPE MESSAGE\ntop file written: \n                      {system_top_file}")
    
    # Create separate itp files for each molecule
    print(f"\nCLEANPIPE MESSAGE\n{len(molecule_names)} itp file written:")
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

    print("\nCLEANPIPE MESSAGE\n\nremoved posres inclusion from: " + s_topfile)

    # Create a more flexible regex pattern to match the inclusion block
    pattern_to_remove = r';\s*Include\s*Position\s*restraint\s*file\s*\n#ifdef\s*POSRES\s*\n#include\s*"posres\.itp"\s*\n#endif\s*'

    # Open the file and read its contents
    with open(s_topfile, 'r') as file:
        file_contents = file.read()

    # Remove the matched block using the regex pattern
    updated_contents = re.sub(pattern_to_remove, '', file_contents, flags=re.MULTILINE)

    # Write the updated contents back to the file
    with open(s_topfile, 'w') as file:
        file.write(updated_contents)




def insert_text_before_directive(s_file_path, s_text_to_insert, s_directive):
    # Read the original content of the file
    with open(s_file_path, 'r') as file:
        lines = file.readlines()

    # Find the line with the specified directive and insert the text before it
    with open(s_file_path, 'w') as file:
        for line in lines:
            if s_directive in line:
                file.write(s_text_to_insert + '\n\n')
            file.write(line)
     

def parse_directive(s_top_file,s_directive_name):
    """
    directives like [ bonds] [ angles ] [ pairs ] etc have several columns, 
    this function returns a list of lists, representing lines and columns nested within each line
   

    """
    list_of_lists = []
    with open(s_top_file, 'r') as file:
        lines = file.readlines()
        inside_constraints = False
        for line in lines:
            if line.strip() == s_directive_name:
                inside_constraints = True
                continue
            if inside_constraints:
                if line.startswith('['):
                    break
                parts_of_line = line.split()
                if len(parts_of_line) >= 2 and parts_of_line[0][0] != ";" and parts_of_line[0][0] != "#":
                    list_of_lists.append(parts_of_line)
                    
    return list_of_lists





def parse_directives_inside_each_and_every_molecule(s_file):
    """
    Parses a GROMACS top/itp file to extract molecule information. The output is something like:

    {
        "SOL": {
            "[ atoms ]": [
                ["1", "SOL", "1", "OW", "0", "0"],
                ["2", "SOL", "1", "HW", "0", "0"],
                ["3", "SOL", "1", "HW", "0", "0"]
            ],
            "[ bonds ]": [
                ["1", "2"],
                ["1", "3"]
            ],
            ...
        },
        "Protein_chain_A": {
            "[ atoms ]": [
                ["1", "MOL", "1", "C", "0.2", "0.0"],
                ["2", "MOL", "1", "O", "0.4", "0.1"]
            ],
            "[ bonds ]": [
                ["1", "2"]
            ],
            ...
        }
    }

    



    The parser only processes sections that belong to a molecule. It assumes that a molecule's data
    starts at a "[ moleculetype ]" directive and ends when either:
      - A new "[ moleculetype ]" directive is encountered, or
      - A system-level directive (like "[ system ]", "[ molecules ]", or "[ intermolecular_interactions ]")
        is encountered.
    
    Parameters:
        s_file (str): The file path to a top or itp
    

    """
    molecules = {}
    current_molecule = None
    current_directive = None
    in_molecule = False

    # Directives that signal the end of molecule-specific information.
    system_directives = {"[ system ]", "[ molecules ]", "[ intermolecular_interactions ]"}
    system_directives_lower = {sd.lower() for sd in system_directives}

    with open(s_file, "r") as file:
        s_file_content = file.read()

    lines = s_file_content.splitlines()
    i = 0
    while i < len(lines):
        # Remove any comment (everything after ';'), any condition statement (everything after '#') and strip whitespace.
        line = lines[i].split(";", 1)[0].split("#", 1)[0].strip()
        #print(line)
        i += 1
        if not line:
            continue


        # Check if the line is a directive (e.g., [ moleculetype ], [ atoms ], etc.)
        if line.startswith("[") and line.endswith("]"):
            directive_lower = line.lower().strip()

            # Start of a new molecule.
            if directive_lower == "[ moleculetype ]":
                in_molecule = True
                current_directive = None  # Reset any previously active sub-directive.

                # Get the molecule header (first non-comment, non-empty line) to obtain the molecule name.
                while i < len(lines):
                    header_line = lines[i].split(";", 1)[0].strip()
                    i += 1
                    if header_line:
                        # The molecule name is taken as the first token.
                        parts = header_line.split()
                        if parts:
                            molecule_name = parts[0]
                            current_molecule = molecule_name
                            molecules[current_molecule] = {}
                        break
                continue

            # If we are within a molecule block and encounter a system-level directive,
            # then we consider the molecule's section to have ended.
            if in_molecule and (directive_lower in system_directives_lower):
                in_molecule = False
                current_molecule = None
                current_directive = None
                continue

            # If we are inside a molecule, then any directive encountered (other than those above)
            # is considered a sub-directive (like "[ atoms ]", "[ bonds ]", etc.).
            if in_molecule:
                current_directive = line  # Store the directive exactly as in the file.
                # Initialize an empty list for this directive if it hasn't been seen yet.
                if current_directive not in molecules[current_molecule]:
                    molecules[current_molecule][current_directive] = []
                continue

            # If not in a molecule, we skip non-relevant directives.
            continue

        # Process non-directive lines: these are data lines that belong to the current sub-directive.
        if in_molecule and current_directive and current_molecule is not None:
            # Split the line into columns based on whitespace.
            row = line.split()
            if row:
                molecules[current_molecule][current_directive].append(row)

    return molecules


def parse_directives_inside_intermolecular_interactions(s_file):
    """
    Parses a GROMACS top/itp file to extract information inside the just directive [ intermolecular_interactions]. 
    
    This is like the function parse_directives_inside_each_and_every_molecule(s_file), but just for the "molecule" [ intermolecular_interactions]
    
    The output is something like:

    {
        "intermolecular_interactions": {
            "[ atoms ]": [
                ["1", "SOL", "1", "OW", "0", "0"],
                ["2", "SOL", "1", "HW", "0", "0"],
                ["3", "SOL", "1", "HW", "0", "0"]
            ],
            "[ bonds ]": [
                ["1", "2"],
                ["1", "3"]
            ],
            ...
        }
    }

    



    The parser only processes whats inside the "[ intermolecular_interactions ]" directive and ends when:
      - A system-level directive (like "[ system ]", "[ molecules ]") is encountered.
    
    Parameters:
        s_file (str): The file path to a top or itp
    

    """
    molecules = {}
    current_molecule = None
    current_directive = None
    in_molecule = False

    # Directives that signal the end of molecule-specific information.
    system_directives = {"[ system ]", "[ molecules ]"}
    system_directives_lower = {sd.lower() for sd in system_directives}

    with open(s_file, "r") as file:
        s_file_content = file.read()

    lines = s_file_content.splitlines()
    i = 0
    while i < len(lines):
        # Remove any comment (everything after ';'), any condition statement (everything after '#') and strip whitespace.
        line = lines[i].split(";", 1)[0].split("#", 1)[0].strip()
        #print(line)
        i += 1
        if not line:
            continue


        # Check if the line is a directive (e.g., [ moleculetype ], [ atoms ], etc.)
        if line.startswith("[") and line.endswith("]"):
            directive_lower = line.lower().strip()

            # Start of [ intermolecular_interactions ].
            if directive_lower == "[ intermolecular_interactions ]":
                in_molecule = True
                current_directive = None  # Reset any previously active sub-directive.

                # define the name that usualy would be the molecule name as "intermolecular_interactions"
                current_molecule = "intermolecular_interactions"
                molecules[current_molecule] = {}
                continue

            # If we are within a molecule block and encounter a system-level directive,
            # then we consider the molecule's section to have ended.
            if in_molecule and (directive_lower in system_directives_lower):
                in_molecule = False
                current_molecule = None
                current_directive = None
                continue

            # If we are inside [ intermolecular_interactions ], then any directive encountered (other than those above)
            # is considered a sub-directive (like "[ atoms ]", "[ bonds ]", etc.).
            if in_molecule:
                current_directive = line  # Store the directive exactly as in the file.
                # Initialize an empty list for this directive if it hasn't been seen yet.
                if current_directive not in molecules[current_molecule]:
                    molecules[current_molecule][current_directive] = []
                continue

            # If not in a molecule, we skip non-relevant directives.
            continue

        # Process non-directive lines: these are data lines that belong to the current sub-directive.
        if in_molecule and current_directive and current_molecule is not None:
            # Split the line into columns based on whitespace.
            row = line.split()
            if row:
                molecules[current_molecule][current_directive].append(row)

    return molecules


import os
import re
from pathlib import Path

# Accept both quotes and angle brackets
INCLUDE_RE = re.compile(r'^\s*#include\s*[<"]([^">]+)[>"]')

def _resolve_include(including_file: Path, include_token: str, search_dirs=()):
    """
    Resolve include_token robustly:
      1) expand ~ and env vars
      2) if absolute -> use as is
      3) if relative -> resolve relative to including file's directory
      4) if not found -> try search_dirs (each joined with relative token)
    Returns (resolved_path: Path | None, tried_paths: list[Path])
    """
    base_dir = including_file.resolve().parent

    expanded = os.path.expandvars(os.path.expanduser(include_token.strip()))
    p = Path(expanded)

    tried = []

    # Candidate 1: absolute or relative-to-including-file
    if p.is_absolute():
        cand = p
    else:
        cand = base_dir / p
    cand = cand.resolve()
    tried.append(cand)
    if cand.exists():
        return cand, tried

    # Candidate 2+: fallback search directories (only meaningful for relative tokens)
    if not p.is_absolute():
        for d in search_dirs:
            d = Path(os.path.expandvars(os.path.expanduser(str(d)))).resolve()
            alt = (d / p).resolve()
            tried.append(alt)
            if alt.exists():
                return alt, tried

    return None, tried


def expand_includes(
    file_path,
    *,
    search_dirs=None,
    visited=None,
    on_missing="keep",   # "keep" or "raise"
    emit_markers=False   # True to add begin/end comments around inserted blocks
):
    """
    Recursively read a file and expand any #include directives.

    Handles:
      - #include "file" and #include <file>
      - absolute + relative include paths
      - ~ and $VARS expansions
      - trailing ';' comments after include directives
      - fallback search directories
      - cycle detection to avoid infinite recursion

    Params:
      search_dirs: iterable of directories to try if relative include isn't found next to including file.
                  Common choices: [os.environ.get("GMXLIB"), "/path/to/ff", os.getcwd()]
      visited: internal set used for cycle detection (you normally don't pass this)
      on_missing: "keep" keeps the original include line; "raise" throws FileNotFoundError
      emit_markers: if True, wraps expanded content with helpful comments.
    """
    if search_dirs is None:
        search_dirs = []
        # Optional: include GMXLIB if set
        gmxl = os.environ.get("GMXLIB")
        if gmxl:
            search_dirs.append(gmxl)

    if visited is None:
        visited = set()

    file_path = Path(file_path).expanduser()
    # Don't force resolve here if it may not exist; but for cycle detection we want a stable key
    file_key = str(file_path.resolve()) if file_path.exists() else str(file_path)

    if file_key in visited:
        # Include loop detected; keep line as-is or raise, but don't recurse forever.
        # Here: keep silently with a marker.
        return f"; NOTE: include loop avoided for {file_path}\n"

    if not file_path.exists():
        msg = f"Top/itp file not found: {file_path}"
        if on_missing == "raise":
            raise FileNotFoundError(msg)
        return f"; NOTE: {msg}\n"

    visited.add(file_key)

    content_parts = []
    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            # Remove trailing GROMACS comments to parse includes robustly,
            # but keep the original line if we don't expand it.
            no_comment = line.split(";", 1)[0].strip()

            m = INCLUDE_RE.match(no_comment)
            if not m:
                content_parts.append(line)
                continue

            include_token = m.group(1).strip()

            resolved, tried = _resolve_include(file_path, include_token, search_dirs=search_dirs)

            if resolved is None:
                tried_str = ", ".join(str(p) for p in tried)
                msg = f"File to be included not found: {include_token}. Tried: {tried_str}"
                if on_missing == "raise":
                    raise FileNotFoundError(msg)
                # Keep original line (your old behavior), but annotate for debugging
                content_parts.append(f"; NOTE: {msg}\n")
                content_parts.append(line)
                continue

            if emit_markers:
                content_parts.append(f"; >>> BEGIN include {include_token} (resolved: {resolved})\n")

            # Recurse
            content_parts.append(
                expand_includes(
                    resolved,
                    search_dirs=search_dirs,
                    visited=visited,
                    on_missing=on_missing,
                    emit_markers=emit_markers,
                )
            )

            if emit_markers:
                content_parts.append(f"; <<< END include {include_token}\n")

    visited.remove(file_key)
    return "".join(content_parts)



def expand_includes_to_temp_file(file_path):
    """
    Expands the given topology file (and all its included files) into one string,
    writes that string to a temporary file (a "virtual file"), and returns the path to that file
    
    Returns:
        A string with the path to the temporary file containing the full, expanded topology.
    """
    # Expand the includes to get the full text
    expanded_text = expand_includes(file_path)
    
    # Create a temporary file to store the expanded top file.
    # Using delete=False so the file remains available via its filename.
    tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.top', delete=False)
    tmp_file.write(expanded_text)
    tmp_file.close()

    print(f"CLEAN PIPE created temp file {tmp_file.name} to store the expanded top file.")
    
    # Return the filename (as a string) that points to the virtual file.
    return tmp_file.name


def deconstruct_top_into_molecules(top_file_path):
    """
    Reads a GROMACS topology file and returns a dictionary mapping each molecule to a path of molecule itp that was generated
    
    Example return value:
    {
    'Protein_chain_A':   'C:\\Users\\Henrique\\AppData\\Local\\Temp\\tmp7dp5w3t9.itp',
     'Support_chain_B':  'C:\\Users\\Henrique\\AppData\\Local\\Temp\\tmpkqp66vn1.itp',
     'SOL':              'C:\\Users\\Henrique\\AppData\\Local\\Temp\\tmp6t421i1v.itp'
     }
    

    The function is smart, as it will expand all the #include statements and extract eache complete 
    [ moleculetype ] blocks – including all lines (comments, blank lines, etc.) – 
    from the first line containing "[ moleculetype ]" until just before a new 
    "[ moleculetype ]" header or a system-level directive header (e.g., 
    "[ system ]", "[ molecules ]", or "[ intermolecular_interactions ]").
    
    The molecule name is determined from the first non-empty, non-comment line 
    immediately following the "[ moleculetype ]" header. Meanwhile, the "[ molecules ]" 
    section is parsed (ignoring blank lines and comment lines) to extract the molecule count.


    example:
    s_top_file = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.top"
    d_mols = cl.deconstruct_top_into_molecules(s_top_file)

    
    """
    # Expand the topology file so that any #include directives are resolved.
    full_top_text = expand_includes(top_file_path)
    lines = full_top_text.splitlines()


    # Regular expression to detect a header line. We check only left-stripped lines.
    header_regex = re.compile(r'^\[\s*(.*?)\s*\]')
    
    # --------- Parse the [ molecules ] section for molecule counts ---------
    molecules_counts = {}
    current_section = None
    for line in lines:
        lstripped = line.lstrip()

        # Check if this line is a header.
        header_match = header_regex.match(lstripped) if lstripped.startswith('[') else None
        if header_match:
            current_section = header_match.group(1).strip().lower()
            continue
        # Only process lines for the molecules section.
        if current_section == "molecules":
            # For count extraction, ignore blank lines and comment lines.
            if not line.strip() or lstripped.startswith(';'):
                continue
            tokens = line.strip().split()
            if len(tokens) >= 2:
                mol_name = tokens[0]
                try:
                    count = int(tokens[1])
                except ValueError:
                    count = tokens[1]  # fallback if not an integer
                molecules_counts[mol_name] = count




    

    # Dictionary to hold complete [ moleculetype ] blocks.
    # Key: molecule name; Value: block content as a single string.
    moleculetype_blocks = {}

    # Regular expression to detect a header line. We check only left-stripped lines.
    header_regex = re.compile(r'^\[\s*(.*?)\s*\]')
    # Headers that terminate a moleculetype block (aside from another moleculetype header)
    termination_headers = {"system", "molecules", "intermolecular_interactions"}

    # Variables to track the current moleculetype block.
    in_moleculetype = False
    current_block = []      # List of lines (preserved exactly as in the file).
    current_mol_name = None # Will be set when the first non-empty, non-comment line is seen.

    # --- Extract complete moleculetype blocks ---
    for line in lines:
        # For header detection we use the left-stripped version of the line.
        lstripped = line.lstrip()
        header_match = header_regex.match(lstripped) if lstripped.startswith('[') else None

        if header_match:
            # We have a header. Extract its content and convert to lowercase.
            header = header_match.group(1).strip().lower()
            # If we are currently in a moleculetype block...
            if in_moleculetype:
                # If the new header indicates the start of a new moleculetype block or
                # is a system-level directive, then we finish the current block.
                if header == "moleculetype" or header in termination_headers:
                    if current_mol_name:
                        # Save the block (all lines collected so far).
                        moleculetype_blocks[current_mol_name] = "\n".join(current_block)
                    # Reset the block state.
                    in_moleculetype = False
                    current_block = []
                    current_mol_name = None
                    # If the header is a new moleculetype header, then start a new block.
                    if header == "moleculetype":
                        in_moleculetype = True
                        current_block.append(line)  # include the header line itself
                    # For a termination header, we do not start a new block.
                    continue
            else:
                # Not currently in a moleculetype block.
                if header == "moleculetype":
                    # Start a new moleculetype block.
                    in_moleculetype = True
                    current_block = [line]  # add the header line
                    current_mol_name = None
                    continue
            # If we're in a moleculetype block and the header is not terminating, we include it.
            if in_moleculetype:
                current_block.append(line)
            # Update current section for molecules count if needed (handled later)
            continue

        # For non-header lines:
        if in_moleculetype:
            # If we haven't yet determined the molecule name, check if this line is
            # not blank and not a comment. (This is only for the purpose of naming the block.)
            if current_mol_name is None:
                if line.strip() and not lstripped.startswith(';'):
                    tokens = line.strip().split()
                    if tokens:
                        current_mol_name = tokens[0]
            # Always add the line to the current block exactly as is.
            current_block.append(line)

    # After processing all lines, if we are still in a moleculetype block, save it.
    if in_moleculetype and current_mol_name:
        moleculetype_blocks[current_mol_name] = "\n".join(current_block)



    # --- Create the final dictionary mapping molecule names to temp_file_path and  count
    molecule_dict = {}
    #these 3 values will be used in just in the first iteration. the goal is to store info of previsous iterations to discover molecule the first_id
    previous_first_id = 0 
    previous_mol_count = 1 
    n_atoms_in_previous_molecule = 1
    for mol_name, count in molecules_counts.items():
        if mol_name in moleculetype_blocks:
            block_content = moleculetype_blocks[mol_name]
            # Write the complete block (including comments and blank lines) to a temporary file.
            tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.itp', delete=False)
            tmp_file.write(block_content)
            tmp_file.close()
            print(f"CLEAN PIPE created temp file {tmp_file.name} to store the itp of the molecule {mol_name}")


            #CONSTRUCT THE DICTIONARY ENTRY
            molecule_dict[mol_name] = tmp_file.name

        else:
            print(f"Warning: Molecule '{mol_name}' declared in [ molecules ] but no [ moleculetype ] block was found.")

    return molecule_dict





def basic_infos_of_molecules(top_file_path):
    """
    Reads top file to obtain total counts of each and all molecules there, and also calculate what will be the fist id.
    
    to do so, it will perform all inclusions, so to add all itp of all molecules to that top, and then parse molecule information

    Example return value:
    {
    'Protein_chain_A':  {'count': 1,   'first_id': 1},
     'Support_chain_B': {'count': 1,   'first_id': 82},
     'SOL':             {'count': 845, 'first_id': 106}
     }
    
    Note: The entire molecule block is written exactly as it appears in the 
    expanded topology file.
    """


    full_top_text = expand_includes(top_file_path)
    lines = full_top_text.splitlines()

    # Regular expression to detect a header line. We check only left-stripped lines.
    header_regex = re.compile(r'^\[\s*(.*?)\s*\]')
    
    # --------- Parse the [ molecules ] section for molecule counts ---------
    molecules_counts = {}
    current_section = None
    for line in lines:
        lstripped = line.lstrip()

        # Check if this line is a header.
        header_match = header_regex.match(lstripped) if lstripped.startswith('[') else None
        if header_match:
            current_section = header_match.group(1).strip().lower()
            continue
        # Only process lines for the molecules section.
        if current_section == "molecules":
            # For count extraction, ignore blank lines and comment lines.
            if not line.strip() or lstripped.startswith(';'):
                continue
            tokens = line.strip().split()
            if len(tokens) >= 2:
                mol_name = tokens[0]
                try:
                    count = int(tokens[1])
                except ValueError:
                    count = tokens[1]  # fallback if not an integer
                molecules_counts[mol_name] = count






    
    # Expand the topology file so that any #include directives are resolved.

    s_top_with_inclusions = expand_includes_to_temp_file(top_file_path)
    ddl_mols_infos = parse_directives_inside_each_and_every_molecule(s_top_with_inclusions)
    

    # --------- Create the final dictionary mapping molecule names to temp_file_path and  count ---------
    molecule_dict = {}
    #these 3 values will be used in just in the first iteration. the goal is to store info of previsous iterations to discover molecule the first_id
    previous_first_id = 0 
    previous_mol_count = 1 
    n_atoms_in_previous_molecule = 1
    for mol_name, directives in molecules_counts.items():

        ll_atoms = ddl_mols_infos[mol_name]['[ atoms ]']# obtain list of lists representing the directive [ atoms ]
        n_atoms_inmolecule = len(ll_atoms) # find the number of atoms
        n_index_last_atom = int(ll_atoms[n_atoms_inmolecule-1][0]) # find the index of the last atom. in gromacs, remember that the index starts at 1
        first_id = previous_first_id + previous_mol_count*n_atoms_in_previous_molecule #calculate the first id. it will use infos of the last iteration

        #CONSTRUCT THE DICTIONARY ENTRY
        molecule_dict[mol_name] = {'count': molecules_counts[mol_name], 'qt_atoms': n_atoms_inmolecule ,'first_id':first_id}

        #save infos that will be used in the next iteration to discover the first_id
        previous_first_id = first_id  
        previous_mol_count =  molecules_counts[mol_name]
        n_atoms_in_previous_molecule = n_atoms_inmolecule



    bricksFileSystem.delete(s_top_with_inclusions)
    return molecule_dict


def get_atom_global_id(s_top_file,s_molecule,n_molecule_instantiation, n_residue, s_atom):

    """
    the input is a top file, a molecule name (eg 'Protein_chain_A'), a resiude id inside that molecule (eg. 4), and an atom name inside that residue (eg CA)
    the function will return the the global id of that atom. this id is the one that will identify the atom insinde the gro file

    example
    
    cl.get_atom_global_id(r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.top",'Protein_chain_A',1,2, 'CA')

    this will return
    5
    """

    s_top_with_inclusions = expand_includes_to_temp_file(s_top_file)

    ddll_parsed_mols = parse_directives_inside_each_and_every_molecule(s_top_with_inclusions)

    #this is the lookup table from wich the local id will be obtained
    ll_atoms = ddll_parsed_mols[s_molecule]['[ atoms ]']

    #this is the borring table that will be filled with the local id
    ll_borring_with_one_line = [[str(n_residue),s_atom]]


    #lookup to find the local id. the key is contructed using the residue number (e.g. 4) and atom name (e.g. CA)
    try:
        ll_filled = lltools.procv_ll(ll_borring_with_one_line, [0,1], ll_atoms, [2,4], [0])
    except:
        print(f'atom not found in [ atoms ] directive of molecule "{s_molecule}"')

    #get the local id what was added in the last column of the borring list of lists
    n_id = int(ll_filled[0][2])


    #get molecule first id, and use it to convert the local id to the global id
    dd_mols = basic_infos_of_molecules(s_top_file)
    n_first_id = dd_mols[s_molecule]['first_id']
    n_qt_atoms = dd_mols[s_molecule]['qt_atoms']
    n_id_global = (n_first_id-1) + (n_molecule_instantiation-1)*n_qt_atoms + n_id
    

    bricksFileSystem.delete(s_top_with_inclusions)
    return n_id_global


def discover_molecule_name_from_global_id(s_top, n_id_global):
    """
    given a gobal id of a atoms, as the ones in the gro file, and a top file,
    this function will return the molecule name that contains that atom
    
    
    
    """





    dd_mols_infos = basic_infos_of_molecules(s_top)
    #this dict contains something like, for example:
    #{ 'Protein_chain_A': {'count': 1, 'qt_atoms': 81, 'first_id': 1},
    #  'Support_chain_B': {'count': 1, 'qt_atoms': 24, 'first_id': 82},
    #  'SOL':             {'count': 845, 'qt_atoms': 3, 'first_id': 106}}

    for key, value in dd_mols_infos.items():
        if n_id_global >= value['first_id'] and n_id_global < (value['first_id'] + value['qt_atoms']*value['count']):
            return key

    # Raise an error if no molecule matches
    raise ValueError(f"No molecule name found for global ID: {n_id_global}")



    return new_ll

def add_lines_at_the_end_of_directive(s_top_file,s_top_file_out, s_directive, ll_lines_to_add, directive_position='first'):
    """
    Inserts a table (list of lists) at the end of a specific directive in a top file.
    
    If more than one occurrence of the directive is found, a warning is printed and the function
    operates on the occurrence specified by the 'directive_position' argument, which can be either 
    'first' or 'last' (default is 'first').

    For example, if the current content of the chosen directive in the top file is:

    [ bonds ]
    1 2 1 1
    2 3 1 1       
               <- ...here is where the lines will be added

    Example usage:
    ll_lines_to_add = [
        ["1", "2", "1", "1"],
        ["2", "3", "1", "1"]
    ]
    s_directive = "[ bonds ]"
    top_file = "oi.top"
    cl.add_lines_at_the_end_of_directive(s_top_file,s_top_file_out, s_directive, ll_lines_to_add, directive_position='first')
    """
    # Validate the directive_position argument.
    if directive_position not in ('first', 'last'):
        raise ValueError("directive_position must be either 'first' or 'last'.")

    # Read the file's content.
    with open(s_top_file, 'r') as file:
        lines = file.readlines()

    # Find all occurrences of the directive.
    directive_indices = [i for i, line in enumerate(lines) if line.strip() == s_directive.strip()]

    if not directive_indices:
        raise ValueError("Directive not found in the file.")

    if len(directive_indices) > 1:
        print(f"Warning: More than one occurrence of directive '{s_directive}' found. Using the {directive_position} occurrence.")

    # Choose the desired occurrence.
    directive_index = directive_indices[0] if directive_position == 'first' else directive_indices[-1]

    # Determine the insertion point: the first subsequent line that starts a new directive (i.e. [ ... ])
    # or the end of the file if no new directive is found.
    insertion_index = len(lines)
    for i in range(directive_index + 1, len(lines)):
        stripped_line = lines[i].strip()
        if stripped_line.startswith('[') and stripped_line.endswith(']'):
            insertion_index = i
            break

    # Convert every cell in each row to a string
    def format_cell(cell):
        if isinstance(cell, float):
            return f"{cell:.4f}"
        return str(cell)
    
    formatted_ll_lines_to_add = [[format_cell(cell) for cell in row] for row in ll_lines_to_add]
    
    # Join each row into a single line (with a newline at the end)
    new_lines = ["   ".join(row) + "\n" for row in formatted_ll_lines_to_add]

    # Insert the new lines.
    updated_lines = lines[:insertion_index] + new_lines + ["\n"] + lines[insertion_index:]

    # Write the updated content back to the file.
    with open(s_top_file_out, 'w') as file:
        file.writelines(updated_lines)


def replace_all_lines_of_directive(s_top_file,s_top_file_out, s_directive, ll_lines_to_add, directive_position='first'):
    """
    Replaces the lines under a specific directive in a top file with a new table (list of lists).
    
    If more than one occurrence of the directive is found, a warning is printed and the function
    operates on the occurrence specified by the 'directive_position' argument, which can be either 
    'first' or 'last' (default is 'first').

    For example, if the current content of the chosen directive in the top file is:

    [ bonds ] 
                <- ...here is where the lines will be added...
    1 2 1 1     <- ...this line will be deleted
    2 3 1 1     <- ...this line will be deleted 
               
    Example usage:
    ll_lines_to_put_there = [
    ['20', '5', '22', '21', '2'],
    ['22', '20', '24', '23', '2'],
    ['30', '24', '32', '31', '2'],
    ['32', '30', '34', '33', '2'],
    ['79', '64', '81', '80', '2'],
    ['1', '5', '20', '22', '2', -47.04540848888724, '5000'],
    ['20', '22', '24', '30', '2', -57.99694175027113, '5000'],
    ['22', '24', '30', '32', '2', -46.481456565006184, '5000'],
    ['52', '54', '60', '62', '2', -46.956430692430864, '5000'],
    ['60', '62', '64', '79', '2', -57.83680231864055, '5000']]

    s_itp_file     = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water_Protein_chain_A.itp"
    s_itp_file_out = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water_Protein_chain_A_edited.top"
    
    cl.replace_all_lines_of_directive(s_itp_file,s_itp_file_out, '[ dihedrals ]', ll_lines_to_put_there, 'last')
    """
    # Validate the directive_position argument.
    if directive_position not in ('first', 'last'):
        raise ValueError("directive_position must be either 'first' or 'last'.")

    # Read the file's content.
    with open(s_top_file, 'r') as file:
        lines = file.readlines()


    # Find all occurrences of the directive.
    directive_indices = [i for i, line in enumerate(lines) if line.strip() == s_directive.strip()]

    if not directive_indices:
        raise ValueError(f"Directive '{s_directive}' not found in the file '{s_top_file}'.")

    if len(directive_indices) > 1:
        print(f"Warning: More than one occurrence of directive '{s_directive}' found. Using the {directive_position} occurrence.")

    # Choose the desired occurrence.
    directive_index = directive_indices[0] if directive_position == 'first' else directive_indices[-1]

    # Determine the end of the directive block.
    block_end_index = len(lines)
    for i in range(directive_index + 1, len(lines)):
        stripped_line = lines[i].strip()
        if stripped_line.startswith('[') and stripped_line.endswith(']'):
            block_end_index = i
            break

    # Convert every cell in each row to a string
    def format_cell(cell):
        if isinstance(cell, float):
            return f"{cell:.4f}"
        return str(cell)
    
    formatted_ll_lines_to_add = [[format_cell(cell) for cell in row] for row in ll_lines_to_add]
    
    # Join each row into a single line (with a newline at the end)
    new_lines = ["   ".join(row) + "\n" for row in formatted_ll_lines_to_add]

    # Replace the lines between the directive header and the block end.
    updated_lines = lines[:directive_index + 1] + new_lines + ["\n"] + lines[block_end_index:]

    # Write the updated content back to the file.
    with open(s_top_file_out, 'w') as file:
        file.writelines(updated_lines)


def put_lines_at_the_proper_place_of_directive(s_file_to_be_edited, s_out_file_name, s_directive, ll_replacement,directive_position='first'):
    """
    The inputs are a list of lists representing a original parsed directive,
    and a list of lists with replacement items that should update the original list.
    
    the original list will be uptaded so that, if the atom ids are present, that line will be replaced. but
    if the atom ids are not present, they will be added at the end



    for example, 
    
    and this is a ll_replacement with new values, that should go into the original ll:
    [['2', '1', '5', '6', '10'],
     ['4', '1', '5', '6', '9', '0.1', '0.2'],
     ['50', '32', '45', '66', '10']]
    
    
    this is a ll_original:
    [['2', '1', '5', '6', '9'], <-this will be replaced
     ['2', '1', '5', '7', '9'],
     ['2', '1', '5', '20', '9'],
     ['3', '1', '5', '6', '9'],
     ['3', '1', '5', '7', '9'],
     ['3', '1', '5', '20', '9'],
     ['4', '1', '5', '6', '9'], <-this will be replaced
     ['4', '1', '5', '6', '9']]  
                                <-there will also be and addition here, for that element that was not found


    EXAMPLE USAGE:
    cl.put_lines_at_the_proper_place_of_directive('protein_in_water.top', 'protein_in_water_new.top', '[ dihedrals ]', ll_gro_diherals_backbone, 'last')

    
    """


    ll_original = parse_directive(s_file_to_be_edited, s_directive)


    # Check if all inputs are lists of lists
    for ll_input in [ll_original, ll_replacement]:
        if not isinstance(ll_input, list):
            raise ValueError(f"Expected input to be a list of lists, but got {type(ll_input).__name__}.")
        if any(not isinstance(row, list) for row in ll_input):
            raise ValueError("Expected input to be a list of lists. Every element in the input list must also be a list")


    #define the number of columns that define the atoms, for each type of directive
    if s_directive in ['[ bonds ]','[ pairs ]','[ constraints ]','[ distance_restraints ]']:
        n_columns_with_atom_ids = 2
    elif s_directive == '[ angles ]':
        n_columns_with_atom_ids = 3
    elif s_directive in ['[ dihedrals ]','[ dihedral_restraints ]']:
        n_columns_with_atom_ids = 4
    elif s_directive == '[ cmap ]':
        n_columns_with_atom_ids = 5
    elif s_directive == '[ position_restraints ]':
        n_columns_with_atom_ids = 1
    else:
        raise ValueError(f"directive {s_directive} not recognized by the function replace_specific_lines_of_directive")
    
    # Create a copy of the original list to avoid in-place modifications
    new_ll = [row.copy() for row in ll_original]

    
    # Convert replacement list-of-lists into a DataFrame (make sure ll2df is defined)
    df = bricksStorage.ll2df(ll_replacement)
    
    # Create a dictionary mapping key (atom aids as strings) to the row index in new_ll
    table_lookup = {
        tuple(str(item) for item in row[:n_columns_with_atom_ids]): idx
        for idx, row in enumerate(new_ll)
    }
    
    # Iterate over each row in the dataframe to update the table
    for _, df_row in df.iterrows():
        # Create key from the atom ids (converted to strings)
        key = tuple(str(x) for x in df_row.iloc[:n_columns_with_atom_ids])
        # Convert the entire dataframe row into a list
        new_row = list(df_row)
        if key in table_lookup:
            new_ll[table_lookup[key]] = new_row
        else:
            # handle the case where the key is not found
            new_ll.append(new_row)


    #insert the updated directive in the top
    replace_all_lines_of_directive(s_file_to_be_edited,s_out_file_name, s_directive, new_ll, directive_position)


def freeze_phi_psi_dihedrals(s_gro_file,s_top_file, restraining_force,s_molename,s_file_to_be_edited, s_out_file_name):
    """

    the goal is to set a dihedral potential in [ dihedrals ] so that the current dihedral angles (according to the first molecule on the top) will remain the same
   
    to acomplish that, for the give molecule, the function will the phi and psi dihedrals, look the gro so to find the real angles,
    and then add those angles as manualy added parameters, togeher with the defined force.
    the funcional is se to be type 2. This is the one that is used for "improper dihedrals", but actually type 2 is just a harmonic potential
    The harmonic potential is the only way to really freeze the dihedral without it being able to rotate because of periodical potentials
    

    the output will be a file where the improper dihedral directive is edited so to add those frozen psi and phi. this file can be either a itp or a top
    depending on where are the dihedral definitions for the chosen molecule

    s_file_to_be_edited. 
    s_out_file_name is 


    s_gro_file             file where the real dihedrals values will be obtained
    s_top_file             file that will inform the ids of the dihedrals of the chosen molecule
    restraining_force      the forceon the harmonic potential. e.g. 5000
    s_molename             the chosen molecule. e.g. 'Protein_chain_A'
    s_file_to_be_edited    the file with information to be edited is. this input exists because sometimes the info is in the top, and sometimes in a itp
    s_out_file_name        just the name of the new file that will be generated. If its the same as s_file_to_be_edited, will override. If its different, the original will be kept
    

    EXAMPLE USAGE:

    s_gro_file     = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.gro"
    s_top_file     = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.top"
    
    restraining_force = 7000
    s_molename = 'Protein_chain_A'
    
    s_file_to_be_edited = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water_Protein_chain_A.itp"
    s_out_file_name = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water_Protein_chain_A_edited.itp"
    
    freeze_phi_psi_dihedrals(s_gro_file,s_top_file, restraining_force,s_molename,s_file_to_be_edited, s_out_file_name)

    """
    
    #parse all the molecules
    s_top_file_with_inclusions = expand_includes_to_temp_file(s_top_file)
    dll_parsed_molecules = parse_directives_inside_each_and_every_molecule(s_top_file_with_inclusions)
    bricksFileSystem.delete(s_top_file_with_inclusions)

    #get dihedralls for the molecule of interest
    ll_dihedrals = dll_parsed_molecules[s_molename]['[ dihedrals ]']
    _ , ll_dihedrals_improper           = lltools.split_ll_diherals_into_proper_and_improper(ll_dihedrals) 

    #get also the atoms, to be able to find each atom name
    ll_atoms     = dll_parsed_molecules[s_molename]['[ atoms ]']

    #extract all dihedrals from gro
    ll_gro_diherals            = bricksGRO.extract_all_dihedrals_from_gro(s_gro_file,s_top_file,s_molename)

    #and select just the ones in the backbone
    ll_gro_diherals_backbone   = lltools.filter_dihedrals_to_keep_only_phi_and_psi(ll_gro_diherals ,ll_atoms)
    
    # organize the list with the backbone dihedrals so we have i,j,k,l,functional,angle,force. all in the correct order
    for row in ll_gro_diherals_backbone:
        angle =row[4]
        row[4]='2'  #will se the functional to 2 (harmonic potentials, suited for improper dihedrals, but also the best option to never twist)
        row.append(angle)# put the angle after the functional
        row.append(str(restraining_force)) #add the restraining force
    
    #add the backbone atoms to the improper dihedrals list. they were in the proper dihedral list and remain there. but now they will also apper in the improper dihedral list
    add_lines_at_the_end_of_directive(s_file_to_be_edited,s_out_file_name, '[ dihedrals ]', ll_gro_diherals_backbone, directive_position='last')

    #update the the dihedral list, but now the atoms of the backbone have improper dihedrals   
    #put_lines_at_the_proper_place_of_directive(s_file_to_be_edited, s_out_file_name, '[ dihedrals ]', ll_gro_diherals_backbone, 'first')
    



def parse_itp_charges(itp_path):
    """
    Parse the [ atoms ] section of a GROMACS .itp and return charges in atom-index order.
    Returns:
        charges: list[float]
        atoms:   list[dict] with keys: nr, atomname, resname, charge

    """
    itp_path = Path(itp_path)
    text = itp_path.read_text(encoding="utf-8", errors="replace").splitlines()

    in_atoms = False
    atoms = []

    for raw in text:
        line = raw.strip()
        if not line:
            continue

        # strip ';' comments
        if ";" in line:
            line = line.split(";", 1)[0].strip()
        if not line:
            continue

        # section headers
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[]").strip().lower()
            in_atoms = (section == "atoms")
            continue

        if not in_atoms:
            continue

        parts = line.split()
        if not parts:
            continue
        if not re.match(r"^\d+$", parts[0]):  # atom index
            continue

        # Typical [ atoms ] format:
        # nr type resnr residue atom cgnr charge mass
        # 0  1    2     3       4    5    6      7
        try:
            nr = int(parts[0])
            resname = parts[3]
            atomname = parts[4]
            charge = float(parts[6])
        except (IndexError, ValueError) as e:
            raise ValueError(f"Failed parsing line in [ atoms ]:\n{raw}") from e

        atoms.append({"nr": nr, "resname": resname, "atomname": atomname, "charge": charge})

    if not atoms:
        raise ValueError(f"No atoms found in [ atoms ] section of {itp_path}")

    atoms.sort(key=lambda d: d["nr"])
    charges = [a["charge"] for a in atoms]
    return charges, atoms




def add_itp_partial_charges_to_bfactor_in_pdb(s_itp_file, pdb_in, pdb_out, field="bfactor", decimals=2):
    """

    - reads charges from itp
    - adds them into preexistent pdb as B-factor
    - prints stats

    example usage
    out_pdb = add_itp_partial_charges_to_bfactor_in_pdb(
        itp_path="CHYO_lipid.itp",
        pdb_in="CHYO.pdb",
        pdb_out="CHYO_charges2.pdb"
    )

    """
    charges, atoms = parse_itp_charges(s_itp_file)

    qsum = sum(charges)
    qmin = min(charges)
    qmax = max(charges)

    print(f"Parsed {len(charges)} charges from: {s_itp_file}")
    print(f"Charge stats: min={qmin:+.4f} e  max={qmax:+.4f} e  sum={qsum:+.4f} e")
    print(f"Writing charges into PDB field: {field}")

    # small preview
    #preview=10
    #print("\nPreview (first few atoms from .itp):")
    #for a in atoms[:preview]:
    #    print(f"  {a['nr']:5d}  {a['resname']:<6s}  {a['atomname']:<6s}  {a['charge']:+.4f}")


    pdb_in = Path(pdb_in)
    pdb_out = Path(pdb_out)

    field = field.lower()
    if field not in ("bfactor", "occupancy"):
        raise ValueError("field must be 'bfactor' or 'occupancy'")

    # PDB fixed columns (0-based slices; end exclusive)
    occ_slice = (54, 60)   # cols 55-60
    bfac_slice = (60, 66)  # cols 61-66

    lines = pdb_in.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_line_indices = [i for i, ln in enumerate(lines) if ln.startswith("ATOM") or ln.startswith("HETATM")]

    if len(atom_line_indices) != len(charges):
        raise ValueError(
            f"Atom count mismatch:\n"
            f"  PDB ATOM/HETATM lines: {len(atom_line_indices)}\n"
            f"  charges provided:       {len(charges)}\n\n"
            f"Fix: make sure pdb_in contains ONLY the atoms that correspond to the .itp, in the same order."
        )

    def fmt(v):
        # PDB occupancy/B-factor are width 6 with typically 2 decimals
        s = f"{v:6.{decimals}f}"
        return s[:6] if len(s) > 6 else s

    for k, idx in enumerate(atom_line_indices):
        ln = lines[idx]
        # ensure length
        if len(ln) < 66:
            ln = ln.ljust(66)

        v = charges[k]
        if field == "occupancy":
            ln = ln[:occ_slice[0]] + fmt(v) + ln[occ_slice[1]:]
        else:
            ln = ln[:bfac_slice[0]] + fmt(v) + ln[bfac_slice[1]:]

        lines[idx] = ln

    pdb_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return pdb_out