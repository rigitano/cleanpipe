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

def check_folder_solvent_box(folder_path):
    # Check if the path is a valid folder
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Error: The folder '{folder_path}' does not exist.")
    
    # Get the folder name
    folder_name = os.path.basename(os.path.normpath(folder_path))
    
    # Create the expected file names
    gro_file = os.path.join(folder_path, folder_name + ".gro")
    top_file = os.path.join(folder_path, folder_name + ".top")
    
    # Check if both .gro and .top files exist with the same base name as the folder
    if os.path.isfile(gro_file) and os.path.isfile(top_file):
        return True
    else:
        raise ValueError(f"Missing either .gro or .top file in the folder '{folder_path}' ")