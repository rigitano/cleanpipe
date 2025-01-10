import os
import subprocess
import sys
import select
import platform
import multiprocessing
import psutil




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

def get_file_location(file_path):
    return os.path.dirname(file_path)

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



def run_and_capture_old(command):
    """
    this funcion will run commands in cmd in a way whats printable in juyter and storable in the output
    """

    print(f"\nCLEANPIPE MESSAGE executing command:\n{command}\n")

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
        print(f"\nCLEANPIPE MESSAGE command failed with return code {process.returncode}")
        if captured_error:
            print(f"\nCLEANPIPE MESSAGE standard error output:\n{captured_error}")
    
    return captured_output + captured_error




def run_and_capture(command):

    """
    Executes a shell command, captures both stdout and stderr in real-time, and waits for it to complete.
    """
    print(f"\nCLEANPIPE MESSAGE executing command:\n{command}\n")

    # Start the subprocess
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1  # Line-buffered mode for real-time capture
    )

    captured_output = ""
    captured_error = ""

    # Reading the output and error streams in real-time
    for line in iter(process.stdout.readline, ''):
        sys.stdout.write(line)
        captured_output += line

    for line in iter(process.stderr.readline, ''):
        sys.stderr.write(line)
        captured_error += line

    # Wait for the process to finish
    process.stdout.close()
    process.stderr.close()
    process.wait()

    if process.returncode != 0:
        print(f"\nCLEANPIPE MESSAGE command failed with return code {process.returncode}")
        if captured_error:
            print(f"\nCLEANPIPE MESSAGE standard error output:\n{captured_error}")

    #return captured_output + captured_error




def create_folder(s_folder_name):
    """
    Check if a folder exists, and create it if it doesn't.

    Parameters:
    folder_name (str): The name of the folder to check/create.
    """
    if not os.path.exists(s_folder_name):
        os.makedirs(s_folder_name)
        print(f"CLEANPIPE MESSAGE: Folder '{s_folder_name}' created.")
    else:
        print(f"CLEANPIPE MESSAGE: Using '{s_folder_name}' previouly created.")

def diagnostics():
    """
    cores
    threads
    cluster
    """
    print(f"Number of processors: {multiprocessing.cpu_count()}")
    print(f"Number of threads: {os.cpu_count()}")
    print(f"Cluster: {platform.node()}")
    # Get the available RAM in GB
    free_ram_gb = psutil.virtual_memory().available / (1024**3)
    print(f"Available RAM: {free_ram_gb:.2f} GB")


def concatenate_files(s_file1, s_file2, s_out_file):
    """
    this have the same result as the bas code:
    cat a.txt b.txt > out.txt


    #example
    cl.concatenate_files('a.txt', 'b.txt', 'out.top')   
    
    """


    with open(s_out_file, 'w') as outfile:
        # Write the content of file1 to the output file
        with open(s_file1, 'r') as infile1:
            outfile.write(infile1.read())
    
    with open(s_out_file, 'a') as outfile:
        # Append the content of file2 to the output file
        with open(s_file2, 'r') as infile2:
            outfile.write(infile2.read())


