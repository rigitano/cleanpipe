import os
import subprocess
import sys
import select
import platform
import multiprocessing
import psutil
import shutil
from pathlib import Path
import threading



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
    file_name, _ = os.path.splitext(file_path) # os.path.splitext() returns a tuple: (filename, extension)
    return os.path.basename(file_name)  # os.path.basename ensures we only get the filename, not the full path

def get_filename_with_extension(file_path):
    file_name, file_extension = os.path.splitext(file_path) # os.path.splitext() returns a tuple: (filename, extension)
    return os.path.basename(file_name) + file_extension  # os.path.basename ensures we only get the filename, not the full path


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

def delete_old(s_filename):
    """
    xxx should be replaced by the new delete function

    example
    delete("posres.itp")
    """
    subprocess.run(f"rm {s_filename}" , shell=True, check=True) 



def delete(target_path):
    """
    Deletes a file or folder safely, working across all operating systems.
    Handles errors like non-existent paths, permission issues, and read-only files.

    # Example usage
    delete_any_path("C:/Users/Example/Desktop/test.txt")  # Windows file
    delete_any_path("/home/user/example.txt")  # Linux/macOS file
    delete_any_path("~/Documents/test_folder")  # Expands to home directory, deletes folder

    """
    try:
        # Convert string to Path object and resolve the absolute path
        path = Path(target_path).expanduser().resolve()

        if not path.exists():
            print(f"CLEAN PIPE Error: {path} does not exist.")
            return

        if path.is_file():
            path.unlink()
            print(f"CLEAN PIPE Deleted file: {path}")
        elif path.is_dir():
            shutil.rmtree(path)  # Recursively delete the folder and all its contents
            print(f"CLEAN PIPE Deleted folder: {path}")
        else:
            print(f"CLEAN PIPE Error: {path} is neither a file nor a folder.")

    except FileNotFoundError:
        print(f"CLEAN PIPE Error: {path} does not exist.")
    except PermissionError:
        print(f"CLEAN PIPE Error: Permission denied to delete {path}. Try running with admin rights.")
    except Exception as e:
        print(f"CLEAN PIPE Unexpected error: {e}")





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




def run_and_capture_old_version_that_works_but_doesnt_stop_code_that_called_it_if_there_was_an_error(command):

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




def run_and_capture(command):
    """
    Executes a shell command while capturing and displaying stdout/stderr in real-time.
    Adds stream identification prefixes to avoid confusion in Jupyter notebooks.
    """
    # Print and flush command message immediately
    print(f"\nCLEANPIPE MESSAGE ### TERMINAL COMMAND ###: {command}", flush=True)
    #print(f"{command}", flush=True)
    
    # Start the process
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )
    
    # Prepare capture variables and synchronization lock
    captured_output = []
    captured_error = []
    output_lock = threading.Lock()
    
    def stream_reader(stream, capture_list, is_stderr=False):
        """Read lines from a stream and handle output with synchronization."""
        while True:
            line = stream.readline()
            if not line:
                break
            with output_lock:
                # Add stream identification prefix
                prefix = "STDERR: " if is_stderr else "STDOUT: "
                full_line = prefix + line
                
                # Write to appropriate sys stream
                if is_stderr:
                    sys.stderr.write(full_line)
                    sys.stderr.flush()
                else:
                    sys.stdout.write(full_line)
                    sys.stdout.flush()
                # Capture the original line
                capture_list.append(line)

    # Create and start threads for stdout/stderr
    stdout_thread = threading.Thread(
        target=stream_reader,
        args=(process.stdout, captured_output, False)
    )
    stderr_thread = threading.Thread(
        target=stream_reader,
        args=(process.stderr, captured_error, True)
    )
    
    stdout_thread.start()
    stderr_thread.start()
    
    # Wait for threads to finish
    stdout_thread.join()
    stderr_thread.join()
    
    # Ensure process completion
    process.wait()
    
    # Convert captured lines to strings
    final_output = ''.join(captured_output)
    final_error = ''.join(captured_error)
    
    # Print completion message
    print(f"CLEANPIPE MESSAGE exit {process.returncode} (ok)\n", flush=True)
    
    # Handle errors
    if process.returncode != 0:
        print(f"CLEANPIPE MESSAGE exit {process.returncode} (fail)\n")
        if final_error:
            print(f"\nCLEANPIPE MESSAGE final standard error output after failiure:\n{final_error}\n")
        raise subprocess.CalledProcessError(
            process.returncode,
            command,
            output=final_output,
            stderr=final_error
        )
    
    return final_output



def run_and_capture_before_adra(command):
    """
    Executes a shell command, capturing both its normal output (stdout) and error output (stderr) in real-time.
    If the command fails (returns a non-zero code), an exception is raised to stop further execution.
    """

    # Tell the user which command is being executed
    print(f"\nCLEANPIPE MESSAGE executing command:\n{command}\n")

    # Start the command as a separate process.
    # - 'shell=True' lets us run the command as if we typed it into the shell.
    # - 'stdout=subprocess.PIPE' and 'stderr=subprocess.PIPE' tell Python to capture the output and errors.
    # - 'text=True' means we work with text (strings) rather than raw bytes.
    # - 'bufsize=1' makes the output be line-buffered, which helps us print output in real-time.
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    # Create empty strings to save the outputs we capture.
    captured_output = ""
    captured_error = ""

    # Read the standard (normal) output line by line as soon as it is available.
    # The iter() function keeps reading until it finds an empty string (which means there's no more output).
    for line in iter(process.stdout.readline, ''):
        sys.stdout.write(line)    # Print the line immediately to the console.
        captured_output += line   # Also add it to our captured_output string.

    # Read the error output line by line in a similar way.
    for line in iter(process.stderr.readline, ''):
        sys.stderr.write(line)    # Print the error line immediately to the console.
        captured_error += line    # Also add it to our captured_error string.

    # Close the output and error streams now that we're done reading from them.
    process.stdout.close()
    process.stderr.close()

    # Wait for the process to finish executing the command.
    process.wait()

    # Check if the command was successful.
    # A return code of 0 usually means success. Any other number indicates an error.
    if process.returncode != 0:
        # Inform the user that the command failed and show the return code.
        print(f"\nCLEANPIPE MESSAGE command failed with return code {process.returncode}")
        if captured_error:
            print(f"\nCLEANPIPE MESSAGE standard error output:\n{captured_error}")
        # Raise an exception so that the caller of this function will not continue running further code.
        # This exception can be caught by the calling code if needed.
        raise subprocess.CalledProcessError(
            process.returncode,
            command,
            output=captured_output,
            stderr=captured_error
        )

    # If everything went well, you might want to return the captured normal output.
    #return captured_output


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
    this have the same result as the bash code:
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


