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
