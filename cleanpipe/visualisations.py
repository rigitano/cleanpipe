import pandas as pd
from io import StringIO
import pandas as pd
import numpy as np
import socket
import subprocess
import tempfile
import os


import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogLocator
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm

import nglview as nv
import mdtraj as md

import seaborn as sns


def visualize_coordinates(s_coord):
    """
    
    example:
    cl.visualize_coordinates("pepticat.pdb")
    """


    view = nv.show_file(s_coord)
    view.clear()

    view.add_representation('cartoon', selection='protein', color='red')
    view.add_representation('ball+stick', selection='protein')

    view.add_representation('ball+stick', selection='not protein', opacity=0.1)

    view.add_representation('licorice', selection='SOL', color='blue', opacity=0.2)
    view.add_representation('licorice', selection='OCT', color='yellow', opacity=0.2)

    view.add_representation('ball+stick', selection='CL', color='yellow', aspectRatio=10)
    view.add_representation('ball+stick', selection='NA', color='green', aspectRatio=10)




    return view



def visualize_trajectory(s_xtc,s_gro):
    """
    
    example:
    cl.visualize_trajectory("pepticat2_in_water/3_NPT/npt.trr","pepticat2_in_water/2_NVT/nvt.gro")
    """
    trajectory  = md.load(s_xtc, top=s_gro)
    view = nv.show_mdtraj(trajectory)
    view.clear()

    view.add_representation('ball+stick', selection='protein')
    view.add_representation('cartoon', selection='protein')

    view.add_representation('ball+stick', selection='CL', color='yellow', aspectRatio=10)
    view.add_representation('ball+stick', selection='NA', color='green', aspectRatio=10)

    view.add_spacefill('not protein', opacity=0.1)
    
    return view



def plot_hbonds(s_ss_file,s_ps_file,s_pp_file,s_subtitle=""):
    """
    the inputs are 3 xvg files that have to be generated calling the tool 'gmx hbonds' 3 times:
    one to calculate solvent-solvent hbonds, other to calculate protein-solvent hbonds, and other to calculate protein-protein hbons

    the last option input is a subtitle

    example:
    cl.plot_hbonds("analysis/hb_ss.xvg","analysis/hb_ps.xvg","analysis/hb_pp.xvg","hello")
    """
   

    # Load data
    def load_data(file):
        with open(file, 'r') as f:
            filtered_lines = ''.join([line for line in f if not line.startswith(('@', '#'))])
        return pd.read_csv(StringIO(filtered_lines), sep=r'\s+', header=None, names=['Time (ps)', 'Number', 'Number2'])

    dfh1 = load_data(s_ss_file)
    dfh2 = load_data(s_ps_file)
    dfh3 = load_data(s_pp_file)

    #convert ps to ns
    dfh1['Time (ns)'] = dfh1['Time (ps)']/1000
    dfh2['Time (ns)'] = dfh2['Time (ps)']/1000
    dfh3['Time (ns)'] = dfh3['Time (ps)']/1000

    # Calculate averages
    avg1 = dfh1['Number'].mean()
    avg2 = dfh2['Number'].mean()
    avg3 = dfh3['Number'].mean()

    # Add a small offset to handle zero values
    eps = 1e-1
    if avg1 == 0:
        dfh1['Number'] += eps
    if avg2 == 0:
        dfh2['Number'] += eps
    if avg3 == 0:
        dfh3['Number'] += eps

    # Plotting data
    plt.figure(figsize=(8, 4))

    plt.plot(dfh1['Time (ns)'], dfh1['Number'], label='Solvent-Solvent', color='navy')
    plt.axhline(avg1, linestyle='--', color='navy', linewidth=1)
    plt.text(max(dfh1['Time (ns)']) * 1.06, avg1+eps, f'{avg1:.1f}', verticalalignment='center', color='navy')

    plt.plot(dfh2['Time (ns)'], dfh2['Number'], label='Protein-Solvent', color='red')
    plt.axhline(avg2, linestyle='--', color='red', linewidth=1)
    plt.text(max(dfh2['Time (ns)']) * 1.06, avg2+eps, f'{avg2:.1f}', verticalalignment='center', color='red')
    
    plt.plot(dfh3['Time (ns)'], dfh3['Number'], label='Protein-Protein', color='saddlebrown')
    plt.axhline(avg3, linestyle='--', color='saddlebrown', linewidth=1)
    plt.text(max(dfh3['Time (ns)']) * 1.06, avg3+eps, f'{avg3:.1f}', verticalalignment='center', color='saddlebrown')



    # Use SymLogScale for y-axis
    plt.yscale('symlog', linthresh=eps)
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(1.0, 10.0) * eps, numticks=10))
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    #plt.grid(True, which='both', linestyle='-', linewidth=0.5)

    # Labels and title
    plt.title('Hydrogen Bonds Over Time\n'+s_subtitle, fontsize=13)
    plt.xlabel('Time (ns)', fontsize=13)
    plt.ylabel('H bond count\n(log-like scale)', fontsize=13)
    plt.ylim(0.1, max(avg1,avg2,avg3)*10)
    
    all_yticks = plt.gca().get_yticks()# Get all existing y-ticks
    filtered_yticks = all_yticks[2:]# Exclude the firsts, because 0.1 is being considered zero, and that can be confusing 
    plt.gca().set_yticks(filtered_yticks)# Set the y-ticks manually, excluding the firsts
    plt.gca().set_yticklabels([f'{int(y):d}' for y in filtered_yticks]) #Add labels with no decimals

    


    # Legend positioning
    plt.legend(bbox_to_anchor=(1.10, 1), loc='upper left')
    plt.tight_layout()
    #plt.grid(True)
    plt.show()

def plot_sasa(s_file,s_subtitle):
    """
    
    example
    cl.plot_sasa('sasa.xvg','hi')
    """

    # load data
    with open(s_file, 'r') as file:
        filtered_lines1 = ''.join([line for line in file if not line.startswith(('@', '#'))])
    data1 = StringIO(filtered_lines1)
    dfh1 = pd.read_csv(data1, sep=r'\s+', header=None, names=['Time (ps)', 'Area'])


    #convert ps to ns
    dfh1['Time (ns)'] = dfh1['Time (ps)']/1000
    
    # Calculate averages
    avg1 = dfh1['Area'].mean()

    # Set the color palette
    sns.set_palette(['#1abc9c'])
    colors = sns.color_palette()
    color1 = colors[0]  # First color in the palette

    
    # Plotting data from df1
    plt.subplots(figsize=(6.5, 3.6))  # Smaller figure size
    plt.plot(dfh1['Time (ns)'], dfh1['Area'], label='Area')
    plt.axhline(avg1, linestyle='--', color=color1)
    if avg1 > 0.1:
        plt.text(max(dfh1['Time (ns)'])*1.06, avg1, f'{avg1:.1f}', verticalalignment='bottom', color=color1)

    plt.title('Solvent-Accessible Surface Area over Time\n' + s_subtitle, fontsize=13)
    plt.xlabel('Time (ns)', fontsize=13)
    plt.ylabel('Area (nm\\S2\\N)', fontsize=13)
    plt.ylim(0, avg1*1.2)
    #plt.yticks([1, 2, 3, 4, 5,6,7,8,9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,22,23,24,25])

    
    #plt.grid(True)
    plt.show()


def plot_ramachandran(s_rama_file,s_subtitle):
    """
    cl.plot_ramachandran("rama.csv","hello")
    """
    
    # Load the data
    data = pd.read_csv(s_rama_file, header=None, names=['phi', 'psi'])
    
    # Define the boundaries for the secondary structure regions
    alpha_region = {'phi': (-180, -50), 'psi': (-60, 45)}
    beta_region = {'phi': (-180, -50), 'psi': (90, 180)}
    left_alpha_region = {'phi': (50, 180), 'psi': (45, 180)}
    
    # Plot the data
    fig, ax = plt.subplots(figsize=(8, 6))  # Smaller figure size
    
    # Background regions
    ax.fill_betweenx(np.linspace(alpha_region['psi'][0], alpha_region['psi'][1], 100), alpha_region['phi'][0], alpha_region['phi'][1], color='red', alpha=0.7, label='Alpha Region')
    ax.fill_betweenx(np.linspace(beta_region['psi'][0], beta_region['psi'][1], 100), beta_region['phi'][0], beta_region['phi'][1], color='yellow', alpha=0.7, label='Beta Region')

    # Scatter plot
    ax.scatter(data['phi'], data['psi'], s=10, color='black', alpha=0.6, label='Residues')
    
    # Axes limits
    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    
    # Labels and title
    ax.set_xlabel('φ', fontsize=13)#Phi
    ax.set_ylabel('ψ', fontsize=13)#Psi
    ax.set_title('Ramachandran Plot\n' + s_subtitle, fontsize=13)
    
    # Grid
    ax.grid(True, linestyle='-', alpha=0.3)
    
    # Legend outside the plot
    ax.legend(fontsize=9, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.8, 0.6])
    
    # Show plot
    plt.show()


def plot_dssp(s_dssp_file,s_subtitle):

    """
    cl.plot_dssp('dssp.dat','hello')
    """

    with open(s_dssp_file, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()
    
    char_to_index = {
        'H': 0,  # alpha-helix
        'B': 1,  # beta-bridge
        'E': 2,  # extended beta-ladder
        'G': 3,  # 3_10-helix
        'I': 4,  # pi-helix
        'P': 5,  # kappa-helix (poly-proline II)
        'S': 6,  # bend
        'T': 7,  # hydrogen-bonded turn
        '=': 8,  # break
        '~': 9   # coil/loop (no structure)
    }
    labels = ['alpha-helix', 
              'beta-bridge', 
              'extended beta-ladder', 
              '3_10-helix', 
              'pi-helix', 
              'kappa-helix', 
              'bend', 
              'hydrogen-bonded turn', 
              'break', 
              'coil/loop (no structure)']


    #colors = ['black','red', 'green', 'blue', 'yellow', 'orange', 'purple', 'brown', 'pink', 'white', 'gray']

    colors = [
    'black',
    'red',        # alpha-helix
    'orange',     # beta-bridge
    'yellow',     # extended beta-ladder
    'darkred',    # 3_10-helix
    'magenta',    # pi-helix
    'pink',       # kappa-helix
    'cyan',       # bend
    'blue',       # hydrogen-bonded turn
    'purple',       # break
    'grey'       # coil/loop (no structure)
]
    
    
    # Strip newline characters and create a 2D list where each sublist is a list of characters from each line
    data = [list(line.strip()) for line in lines]
    
    # Transpose the data to switch rows and columns
    data_transposed = list(zip(*data))
    
    # Map characters to numbers using the provided dictionary
    mapped_data = [[char_to_index.get(char, -1) for char in row] for row in data_transposed]
    
    # Convert the mapped data into a NumPy array
    data_indices = np.array(mapped_data)
    
    # Add a row of -1 at the beginning and a column of -1 at the beginning of each row, shifting the data, so now indexes of real data start at 1 instead of zero
    data_indices = np.pad(data_indices, ((1, 0), (1, 0)), mode='constant', constant_values=-1)
    #print(data_indices)
    
    # Define colors for each type

    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(-1.5, 10, 1)  # This sets bounds at midpoints between integers from -1 to 9
    
    # Create a BoundaryNorm which will use the specified bounds
    norm = BoundaryNorm(bounds, cmap.N)
    
    # Read the data from the file
    
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(12, 3))
                                   
    #heatmap = ax.imshow(data_indices, aspect='auto', cmap=cmap, norm=norm)
    # Get the dimensions of the matrix
    num_rows, num_cols = data_indices.shape
    heatmap = ax.imshow(data_indices[1:,1:], aspect='auto', cmap=cmap, interpolation='nearest',norm=norm, extent=[0.5, num_cols-0.5, num_rows-0.5, 0.5])
    
    
    # Create colorbar with labels
    #colorbar = plt.colorbar(heatmap, ticks=np.arange(len(colors)))
    #colorbar.set_ticklabels([
    #    'alpha-helix', 'isolated beta-bridge', 'extended strand in beta-ladder',
    #    '3_10-helix', 'pi-helix', 'kappa-helix', 'bend', 'hydrogen-bonded turn', 'break', 'loop'
    #])
    
    

    legend_handles = [mpatches.Patch(color=colors[i+1], label=labels[i]) for i in range(len(labels))] # be carefull, there are more colors than labels, because there is black for -1, and that shouldnt have any label
    legend = plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Set axis labels
    ax.set_xlabel('Time (ns)', fontsize=13)
    ax.set_ylabel('Residue Position', fontsize=13)
    ax.set_title('Secondary Structure Over Time (DSSP Algorithm)\n' + s_subtitle, fontsize=13)
    plt.grid(False)
    
    # Optionally save the mapped data to a file
    #np.savetxt('transposed_data_indices.dat', data_indices, fmt='%d', delimiter=',', header='Transposed Data Indices')

    # Show the plot
    plt.show()


def open_vmd_with_socket():
    """"
    xxx this should be used in the app

    """

    # Define the Tcl script as a string
    tcl_script = """
    proc start_server {port} {
        set server [socket -server handle_connection $port]
        puts "Server started on port $port"
        return $server
    }

    proc handle_connection {sock addr port} {
        puts "Connection from $addr:$port"
        fconfigure $sock -buffering line
        while {[gets $sock line] >= 0} {
            puts "Received command: $line"
            catch {eval $line} result
            puts $sock $result
            flush $sock
        }
        close $sock
    }

    start_server 5555
    """

    # Create a temporary file for the script
    with tempfile.NamedTemporaryFile(delete=False, suffix=".tcl") as temp_script:
        temp_script.write(tcl_script.encode('utf-8'))
        temp_script_path = temp_script.name

    try:
        # Call VMD with the temporary script
        subprocess.run(["C:\\Program Files\\VMD\\vmd", "-e", temp_script_path], check=True)
    finally:
        # Ensure the temporary file is deleted
        if os.path.exists(temp_script_path):
            os.remove(temp_script_path)






def open_vmd_with_socket2():
    """
    Launches VMD with a Tcl script for a socket server without blocking the Jupyter cell.
    """

    # Define the Tcl script as a string
    tcl_script = """
    proc start_server {port} {
        set server [socket -server handle_connection $port]
        puts "Server started on port $port"
        return $server
    }

    proc handle_connection {sock addr port} {
        puts "Connection from $addr:$port"
        fconfigure $sock -buffering line
        while {[gets $sock line] >= 0} {
            puts "Received command: $line"
            catch {eval $line} result
            puts $sock $result
            flush $sock
        }
        close $sock
    }

    start_server 5555
    """

    # Create a temporary file for the script
    with tempfile.NamedTemporaryFile(delete=False, suffix=".tcl") as temp_script:
        temp_script.write(tcl_script.encode('utf-8'))
        temp_script_path = temp_script.name

    try:
        # Call VMD with the temporary script in a separate process
        subprocess.Popen(
            ["C:\\Program Files\\VMD\\vmd", "-e", temp_script_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            close_fds=True
        )
    finally:
        # Ensure the temporary file is deleted
        if os.path.exists(temp_script_path):
            #os.remove(temp_script_path)
            print("ok")



def send_command_to_vmd(s_command):
    """"
    xxx this should be used in the app


    ex: 
    cl.send_command_to_vmd("graphics top cylinder {0 0 0} {10 10 10} radius 0.1") # this command creates a cilinder: 
    """

    print("CLEAN PIPE command sent to vmd:\n" + s_command)

    send_command_to_vmd
    with socket.create_connection(("localhost", 5555)) as sock:
        command = s_command 
        sock.sendall(command.encode('utf-8') + b'\n')
        response = sock.recv(1024)
        print("Response:", response.decode('utf-8'))

