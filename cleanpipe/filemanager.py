import os

def check_file(file_path,v_alowedExtentions):
    '''
    given the name of the file, and the alowed extentions, the function will check if the file is valid
    '''

    if os.path.isfile(file_path):
        # Get the file extension
        _, file_extension = os.path.splitext(file_path)
        # Check if the file has the correct extension (.pdb or .gro)
        if file_extension in v_alowedExtentions:
            return "ok"
        else:
            raise ValueError(f"Error: Invalid file extension '{file_path}'. Only .pdb or .gro are allowed.")
    else:
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")
    
def get_filename_without_extension(file_path):
    # os.path.splitext() returns a tuple: (filename, extension)
    file_name, _ = os.path.splitext(file_path)
    return os.path.basename(file_name)  # os.path.basename ensures we only get the filename, not the full path

def check_folder(folder_path):
    # Check if the path is a valid folder
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Error: The folder '{folder_path}' does not exist.")
    else:
        return True
    

def get_gro_and_top(folder_path):
    gro_files = []
    top_files = []

    # Loop through the files in the given folder
    for file_name in os.listdir(folder_path):
        # Check if the file has .gro extension
        if file_name.endswith(".gro"):
            gro_files.append(file_name)
        # Check if the file has .top extension
        elif file_name.endswith(".top"):
            top_files.append(file_name)

    # Check if there's exactly one .gro and one .top file
    if len(gro_files) != 1:
        raise ValueError(f"Expected exactly one .gro file, found {len(gro_files)}")
    
    if len(top_files) != 1:
        raise ValueError(f"Expected exactly one .top file, found {len(top_files)}")
    
    return gro_files[0], top_files[0]

