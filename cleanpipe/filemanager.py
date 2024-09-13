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
    

def get_all_files_with_certain_extention(s_folder_path,s_extention):
    """
    search a fiven folder to return a list of the file names with a given extention
    """


    l_files = []

    # Loop through the files in the given folder
    for s_file_name in os.listdir(s_folder_path):
        # Check if the file has .gro extension
        if s_file_name.endswith(s_extention):
            l_files.append(s_file_name)
    return l_files


def get_single_gro(s_folder_path):
    """
    get a gro file in a folder. the folder must contain only one gro
    """

    l_files = get_all_files_with_certain_extention(s_folder_path,".gro")

    # Check if there's exactly one .gro and one .top file
    if len(l_files) != 1:
        raise ValueError(f"Expected exactly one .gro file, found {len(l_files)}")
    else:
        return l_files[0]




def get_single_top(s_folder_path):
    """
    get a top file in a folder. the folder must contain only one top
    """

    l_files = get_all_files_with_certain_extention(s_folder_path,".top")

    # Check if there's exactly one .gro and one .top file
    if len(l_files) != 1:
        raise ValueError(f"Expected exactly one .top file, found {len(l_files)}")
    else:
        return l_files[0]


def get_all_itps(s_folder_path):


    """
    get all the itp files in a given folder
    """

    l_files = get_all_files_with_certain_extention(s_folder_path,".itp")

    return l_files