import os
import subprocess
import sys
import select


def check_extention(file_path,v_alowedExtentions):
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

def delete(s_filename):
    """
    example
    delete("posres.itp")
    """
    subprocess.run(f"rm {s_filename}" , shell=True, check=True) 



def run_and_capture(command):
    """
    this funcion will run commands in cmd in a way whats printable in juyter and storable in the output
    """

    print(f"\nCLEANPIPE: executing command {command}\n")

    captured_output = ""
    captured_error = ""

    # Start the subprocess
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # Monitor the stdout and stderr streams
    while True:
        # Use select to check for available output
        readable, _, _ = select.select([process.stdout, process.stderr], [], [])

        for stream in readable:
            # Read a single line from the stream
            line = stream.readline()
            if line == "" and process.poll() is not None:
                # End of stream and process has exited
                break
            if line:
                if stream == process.stdout:
                    sys.stdout.write(line)  # Print stdout to notebook in real-time
                    captured_output += line  # Store the stdout line
                elif stream == process.stderr:
                    sys.stderr.write(line)  # Print stderr to notebook in real-time
                    captured_error += line  # Store the stderr line

        # Break out of the loop if the process is finished
        if process.poll() is not None:
            break

    # Wait for the process to finish and get the return code
    process.wait()

    if process.returncode != 0:
        print(f"CLEAN PIPE: Command failed with return code {process.returncode}")
        if captured_error:
            print(f"CLEAN PIPE: Standard Error Output:\n{captured_error}")
    
    return captured_output + captured_error