import pandas as pd
from io import StringIO
import pandas as pd
import numpy as np
import socket
import subprocess
import tempfile
import os
import platform


from cleanpipe import lltools
from cleanpipe import algelin
from cleanpipe import bricksTOP
from cleanpipe import bricksGRO
from cleanpipe import bricksFileSystem



import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogLocator
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm

import nglview as nv
import MDAnalysis as mda
import mdtraj as md

import seaborn as sns


def view_coord(s_coord):
    """
    
    example:
    cl.view_coord("pepticat.pdb")
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



def view_traj(s_xtc,s_gro):
    """
    
    example:
    cl.view_traj("pepticat2_in_water/3_NPT/npt.trr","pepticat2_in_water/2_NVT/nvt.gro")
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
        system = platform.system()

        if system == "Windows":
            # Windows: direct path to VMD executable
            vmd_command = [
                "C:\\Program Files\\VMD\\vmd",
                "-e", temp_script_path
            ]
            subprocess.Popen(
                vmd_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                close_fds=True
            )

        elif system == "Linux":
            # Linux: run bash as login shell so ~/.bashrc gets sourced
            bash_command = f"source /etc/profile.d/modules.sh && module load vmd && vmd -e {temp_script_path}"

            #env = os.environ.copy()
            #env["DISPLAY"] = ":1"

            #print("Trying to launch VMD with:")
            #print("Command:", bash_command)
            #print("Environment DISPLAY:", env["DISPLAY"])

            
            process = subprocess.Popen(
                ["/bin/bash", "-l", "-c", bash_command],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                close_fds=True
                
            )
            #stdout, stderr = process.communicate(timeout=10)
            #print("STDOUT:", stdout.decode())
            #print("STDERR:", stderr.decode())


        else:
            raise OSError(f"vmd is not in 'C:\\Program Files\\VMD\\vmd' (for windows), nor callable using 'module load vmd' (for linux). This is your system: {system}")
    finally:
        # Ensure the temporary file is deleted
        if os.path.exists(temp_script_path):
            #os.remove(temp_script_path)
            print(f"DISPLAY: {os.environ.get('DISPLAY')}")

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



def see_interactions(s_top,s_gro,s_mol_name):

    """
    based on the top, this funcion will draw the topology.

    each parameter of each directive will be ploted as a new molecule
    parameters that as zero will be ploted as white




    example

    s_top = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/pepticat9_truss_in_water.top"
    s_gro = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/3_NPT/npt.gro"

    cl.see_interactions(s_top,s_gro,'Protein_chain_A')
    
    """





    print("CLEAN PIPE excecuting function see_topology_of_one_molecule")





    
    
    what_to_plot = ["[ bonds ]","[ angles ]","[ dihedrals ]","[ cmap ]","[ pairs ]","[ exclusions ]","[ position_restraints ]"]
    #what_to_plot = ["[ bonds ]","[ dihedrals ]","[ cmap ]","[ position_restraints ]"]
    #what_to_plot = ["[ dihedrals ]","[ cmap ]","[ pairs ]","[ exclusions ]","[ position_restraints ]"]
    #what_to_plot = ["[ bonds ]","[ angles ]","[ dihedrals ]","[ cmap ]","[ pairs ]","[ exclusions ]"]
    #what_to_plot = ["[ position_restraints ]"]
    #what_to_plot = ["[ bonds ]"]
    
    #l_tcl_commads.append("graphics top delete all")

    #here are usefull advanged setups, of transoformations for a give molecule. I fixed the issue in a simpler way, using the command "display resetview"
    #this obtains a transformatin matrix:
    #"molinfo top get center_matrix"
    #this changes the transformation matrices
    #l_tcl_commads.append("molinfo [molinfo top] set {center_matrix rotate_matrix scale_matrix global_matrix} {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}")

    print("CLEAN PIPE gathering topology data")




    #a function to help to deal with items not found in the lookup, probably because of the atom order
    def split_ll_into_found_and_not_found(ll_input):
        #split based on the last element last element equals '2'
        ll_found = [sublist for sublist in ll_input if sublist[-1] != 'Unknown']
        # here we exclude the last column
        ll_not_found = [sublist[:(len(sublist)-1)] for sublist in ll_input if sublist[-1] == 'Unknown']
        return ll_found, ll_not_found


    
    #inser all inclusions of the top file, so that all forcefield information will be there, at the same place
    s_top_with_inclusions = bricksTOP.expand_includes_to_temp_file(s_top)






    # get parameter infos from the top file. please notice that those infos were originally in he ffbonded.itp file, but we included that in the top
    ll_bondtypes     = lltools.clean_comments_out(bricksTOP.parse_directive(s_top_with_inclusions, "[ bondtypes ]"))
    ll_angletypes    = lltools.clean_comments_out(bricksTOP.parse_directive(s_top_with_inclusions, "[ angletypes ]"))
    ll_dihedraltypes = lltools.clean_comments_out(bricksTOP.parse_directive(s_top_with_inclusions, "[ dihedraltypes ]"))


    #now parse each and every molecule in the top
    dd_parsed_mols = bricksTOP.parse_directives_inside_each_and_every_molecule(s_top_with_inclusions)

    #if there are [ intermolecular_interactions ], add it as if it were a molecule. just dont forget that the ids are global there
    dd_parsed_intermolecular = bricksTOP.parse_directives_inside_intermolecular_interactions(s_top_with_inclusions)
    if dd_parsed_intermolecular != {}:
        dd_parsed_mols['intermolecular_interactions'] = dd_parsed_intermolecular['intermolecular_interactions']


    #understand how many molecules there are in the top, and the number of atoms in each
    dd_mols_infos = bricksTOP.basic_infos_of_molecules(s_top_with_inclusions)
    #this dict contains something like, for example:
    #{ 'Protein_chain_A': {'count': 1, 'qt_atoms': 81, 'first_id': 1},
    #  'Support_chain_B': {'count': 1, 'qt_atoms': 24, 'first_id': 82},
    #  'SOL':             {'count': 845, 'qt_atoms': 3, 'first_id': 106}}



    #obtain infos for the molecule of interest
    n_first_id       = dd_mols_infos.get(s_mol_name, {}).get('first_id', 1)  
    n_atoms_in_mol   = dd_mols_infos.get(s_mol_name, {}).get('qt_atoms', 0)  
    n_molecules      = dd_mols_infos.get(s_mol_name, {}).get('count', 1)     

    
    #obtain parsed directives of a certain molecule. Im using gte because some molecule might not have a certain directive. this is not a problem.

    ll_atoms      = dd_parsed_mols.get(s_mol_name, {}).get('[ atoms ]', [])
    ll_bonds      = dd_parsed_mols.get(s_mol_name, {}).get('[ bonds ]', [])
    ll_angles     = dd_parsed_mols.get(s_mol_name, {}).get('[ angles ]', [])
    ll_dihedrals  = dd_parsed_mols.get(s_mol_name, {}).get('[ dihedrals ]', [])
    ll_cmap       = dd_parsed_mols.get(s_mol_name, {}).get('[ cmap ]', [])
    ll_pairs      = dd_parsed_mols.get(s_mol_name, {}).get('[ pairs ]', [])
    ll_exclusions = dd_parsed_mols.get(s_mol_name, {}).get('[ exclusions ]', [])
    ll_posres     = dd_parsed_mols.get(s_mol_name, {}).get('[ position_restraints ]', [])




    


    #now obtain the gro, so to get the coordinates as needed
    ld_coordinates = bricksGRO.parse_gro(s_gro)
    df_coordinates = pd.DataFrame(ld_coordinates)
    df_coordinates.set_index('id', inplace=True)# Set 'id' as the index for fast lookups
    #print(df_coordinates)

    
    #########   add atom names to the table, looking up the [ atoms ] directive ########
    
    if s_mol_name != 'intermolecular_interactions': #within a certain molecule, we just lookup the [ atoms ] directive
        
        # BONDS
        ll_bonds_named = lltools.procv_ll(ll_bonds,[0],ll_atoms,[0],[1])
        ll_bonds_named = lltools.procv_ll(ll_bonds_named,[1],ll_atoms,[0],[1])


        #ANGLES
        ll_angles_named = lltools.procv_ll(ll_angles,[0],ll_atoms,[0],[1])
        ll_angles_named = lltools.procv_ll(ll_angles_named,[1],ll_atoms,[0],[1])
        ll_angles_named = lltools.procv_ll(ll_angles_named,[2],ll_atoms,[0],[1])

        #DIHEDRALS
        ll_dihedrals_named = lltools.procv_ll(ll_dihedrals,[0],ll_atoms,[0],[1])
        ll_dihedrals_named = lltools.procv_ll(ll_dihedrals_named,[1],ll_atoms,[0],[1])
        ll_dihedrals_named = lltools.procv_ll(ll_dihedrals_named,[2],ll_atoms,[0],[1])
        ll_dihedrals_named = lltools.procv_ll(ll_dihedrals_named,[3],ll_atoms,[0],[1])
    
    elif s_mol_name == 'intermolecular_interactions': 
        #it this case its a bit harder, because there's a global id, and no [ atoms ] directive, 
        #so for each global id it will be necessary to find the correct molecule to get its [ atoms ] and calculate the intermolegular directive from the global ones

        #BONDS
        ll_bonds_named = [] #
        for line in ll_bonds: #

            #get the ids in the first couple of columns, but in this case they are global
            global_id1 = int(line[0])
            global_id2 = int(line[1])

            #this function will discover the molecule name using the global ids. this will be usefull to get the correct [ atoms ] directive, where the name is
            current_s_mol_name1 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id1)
            current_s_mol_name2 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id2)
            
            #obtain [ atoms ] directive of the apropriate molecule. the names are there
            current_ll_atoms1      = dd_parsed_mols.get(current_s_mol_name1, {}).get('[ atoms ]', [])
            current_ll_atoms2      = dd_parsed_mols.get(current_s_mol_name2, {}).get('[ atoms ]', [])
            
            #obtain basif infos for the current molecules, this will be usefull do calculate the intermolar ids from the global ones
            n_first_id1       = dd_mols_infos.get(current_s_mol_name1, {}).get('first_id', 1)  
            n_first_id2       = dd_mols_infos.get(current_s_mol_name2, {}).get('first_id', 1)  

            #calculate the intramolecular ids instead of the global ones. make them a str too, because thats how there were stored in the lookup table
            intramolecular_id1 = str(global_id1 - (n_first_id1-1))
            intramolecular_id2 = str(global_id2 - (n_first_id2-1))
            
            #lookput, but I have to create a borring table with just one line, and then I have to extract the first line and last colum of the result to get the string with the atom name
            found_name1 = lltools.procv_ll([[intramolecular_id1]],[0],current_ll_atoms1,[0],[1])[0][1] 
            found_name2 = lltools.procv_ll([[intramolecular_id2]],[0],current_ll_atoms2,[0],[1])[0][1] 

            #concatenate all the information in a single list, and append it to the table
            ll_bonds_named.append(line + [found_name1] + [found_name2]) #


        #ANGLES
        ll_angles_named = [] #
        for line in ll_angles: #

            #get the ids in the first couple of columns, but in this case they are global
            global_id1 = int(line[0])
            global_id2 = int(line[1])
            global_id3 = int(line[2])

            #this function will discover the molecule name using the global ids. this will be usefull to get the correct [ atoms ] directive, where the name is
            current_s_mol_name1 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id1)
            current_s_mol_name2 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id2)
            current_s_mol_name3 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id3)
            
            #obtain [ atoms ] directive of the apropriate molecule. the names are there
            current_ll_atoms1      = dd_parsed_mols.get(current_s_mol_name1, {}).get('[ atoms ]', [])
            current_ll_atoms2      = dd_parsed_mols.get(current_s_mol_name2, {}).get('[ atoms ]', [])
            current_ll_atoms3      = dd_parsed_mols.get(current_s_mol_name3, {}).get('[ atoms ]', [])
            
            #obtain basif infos for the current molecules, this will be usefull do calculate the intermolar ids from the global ones
            n_first_id1       = dd_mols_infos.get(current_s_mol_name1, {}).get('first_id', 1)  
            n_first_id2       = dd_mols_infos.get(current_s_mol_name2, {}).get('first_id', 1)  
            n_first_id3       = dd_mols_infos.get(current_s_mol_name3, {}).get('first_id', 1)  

            #calculate the intramolecular ids instead of the global ones. make them a str too, because thats how there were stored in the lookup table
            intramolecular_id1 = str(global_id1 - (n_first_id1-1))
            intramolecular_id2 = str(global_id2 - (n_first_id2-1))
            intramolecular_id3 = str(global_id3 - (n_first_id3-1))
            
            #lookput, but I have to create a borring table with just one line, and then I have to extract the first line and last colum of the result to get the string with the atom name
            found_name1 = lltools.procv_ll([[intramolecular_id1]],[0],current_ll_atoms1,[0],[1])[0][1] 
            found_name2 = lltools.procv_ll([[intramolecular_id2]],[0],current_ll_atoms2,[0],[1])[0][1] 
            found_name3 = lltools.procv_ll([[intramolecular_id3]],[0],current_ll_atoms3,[0],[1])[0][1] 

            #concatenate all the information in a single list, and append it to the table
            ll_angles_named.append(line + [found_name1] + [found_name2] + [found_name3]) #

        
        #DIHEDRALS
        ll_dihedrals_named = [] #
        for line in ll_dihedrals: #

            #get the ids in the first couple of columns, but in this case they are global
            global_id1 = int(line[0])
            global_id2 = int(line[1])
            global_id3 = int(line[2])
            global_id4 = int(line[4])

            #this function will discover the molecule name using the global ids. this will be usefull to get the correct [ atoms ] directive, where the name is
            current_s_mol_name1 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id1)
            current_s_mol_name2 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id2)
            current_s_mol_name3 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id3)
            current_s_mol_name4 = bricksTOP.discover_molecule_name_from_global_id(s_top_with_inclusions, global_id4)
            
            #obtain [ atoms ] directive of the apropriate molecule. the names are there
            current_ll_atoms1      = dd_parsed_mols.get(current_s_mol_name1, {}).get('[ atoms ]', [])
            current_ll_atoms2      = dd_parsed_mols.get(current_s_mol_name2, {}).get('[ atoms ]', [])
            current_ll_atoms3      = dd_parsed_mols.get(current_s_mol_name3, {}).get('[ atoms ]', [])
            current_ll_atoms4      = dd_parsed_mols.get(current_s_mol_name4, {}).get('[ atoms ]', [])
            
            #obtain basif infos for the current molecules, this will be usefull do calculate the intermolar ids from the global ones
            n_first_id1       = dd_mols_infos.get(current_s_mol_name1, {}).get('first_id', 1)  
            n_first_id2       = dd_mols_infos.get(current_s_mol_name2, {}).get('first_id', 1)  
            n_first_id3       = dd_mols_infos.get(current_s_mol_name3, {}).get('first_id', 1)  
            n_first_id4       = dd_mols_infos.get(current_s_mol_name4, {}).get('first_id', 1)  

            #calculate the intramolecular ids instead of the global ones. make them a str too, because thats how there were stored in the lookup table
            intramolecular_id1 = str(global_id1 - (n_first_id1-1))
            intramolecular_id2 = str(global_id2 - (n_first_id2-1))
            intramolecular_id3 = str(global_id3 - (n_first_id3-1))
            intramolecular_id4 = str(global_id4 - (n_first_id4-1))
            
            #lookput, but I have to create a borring table with just one line, and then I have to extract the first line and last colum of the result to get the string with the atom name
            found_name1 = lltools.procv_ll([[intramolecular_id1]],[0],current_ll_atoms1,[0],[1])[0][1] 
            found_name2 = lltools.procv_ll([[intramolecular_id2]],[0],current_ll_atoms2,[0],[1])[0][1] 
            found_name3 = lltools.procv_ll([[intramolecular_id3]],[0],current_ll_atoms3,[0],[1])[0][1] 
            found_name4 = lltools.procv_ll([[intramolecular_id4]],[0],current_ll_atoms4,[0],[1])[0][1] 

            #concatenate all the information in a single list, and append it to the table
            ll_dihedrals_named.append(line + [found_name1] + [found_name2] + [found_name3] + [found_name4]) #



    else:
        raise ValueError(f"Molecule name not found: {s_mol_name}")
    


    ######## add parameters values to the table, looking up the info that came from the ffbonded file ########

    #BONDS
    ll_bonds_filled = lltools.procv_ll(ll_bonds_named,[set([3,4]),2], ll_bondtypes,[set([0,1]),2], list(range(3,len(ll_bondtypes[0]))))

    #ANGLES
    ll_angles_filled = lltools.procv_ll(ll_angles_named,[set([4,6]),5,3], ll_angletypes,[set([0,2]),1,3], list(range(4,len(ll_angletypes[0]))))

    #DIHEDRALS, they are more complex because of the need of spliting into proper and improper and the possiblility of atom reordering
    
    #split into proper and improper
    ll_dihedraltypes_proper, ll_dihedraltypes_improper       = lltools.split_ll_diherals_into_proper_and_improper(ll_dihedraltypes)
    ll_dihedrals_proper, ll_dihedrals_improper               = lltools.split_ll_diherals_into_proper_and_improper(ll_dihedrals)
    ll_dihedrals_filled_proper, ll_dihedrals_filled_improper = lltools.split_ll_diherals_into_proper_and_improper(ll_dihedrals_named)
        
    # add paramenters for the given atomtypes for PROPER AND IMPROPER
    ll_dihedrals_filled_proper   = lltools.procv_ll(ll_dihedrals_filled_proper,[5,6,7,8,4], ll_dihedraltypes_proper,[0,1,2,3,4], list(range(5,len(ll_dihedraltypes_proper[0]))))# add paramenters for the given atomtypes
    ll_dihedrals_filled_improper = lltools.procv_ll(ll_dihedrals_filled_improper,[5,6,7,8,4], ll_dihedraltypes_improper,[0,1,2,3,4], list(range(5,len(ll_dihedraltypes_improper[0]))))# add paramenters for the given atomtypes
        
    #check the reverse order for cases not found
    
    ll_dihedrals_filled_proper, ll_not_found = split_ll_into_found_and_not_found(ll_dihedrals_filled_proper)
    ll_second_try = lltools.procv_ll(ll_not_found,[8,7,6,5,4], ll_dihedraltypes_proper,[0,1,2,3,4], list(range(5,len(ll_dihedraltypes_proper[0]))))# add paramenters for the given atomtypes
    ll_dihedrals_filled_proper = ll_dihedrals_filled_proper + ll_second_try
    
    ll_dihedrals_filled_improper, ll_not_found = split_ll_into_found_and_not_found(ll_dihedrals_filled_improper)
    
    ll_second_try = lltools.procv_ll(ll_not_found,[8,7,6,5,4], ll_dihedraltypes_improper,[0,1,2,3,4], list(range(5,len(ll_dihedraltypes_improper[0]))))# add paramenters for the given atomtypes
    ll_dihedrals_filled_improper = ll_dihedrals_filled_improper + ll_second_try



    


    ######## rewrite filled list if there are manually defined parameters in the directives  #######
    #for now [ bonds ], [ angles ], [ dihedrals ]
    print("CLEAN PIPE checking for used defined parameters")


    # [ bonds ]
    for i in range(0, len(ll_bonds)): #
        line_unfilled = ll_bonds[i] #
        line_filled   = ll_bonds_filled[i] #

        n_len_line_u = len(line_unfilled)
        n_len_line_f  = len(line_filled)

        #in [ bonds ] Im looking for this type of filled columns in the current line
        #[ bonds ]
        #;  ai    aj funct            c0            c1            c2            c3
        #   1     2     1           0.1234       0.4321
        n_ids_and_functional = 3 # bonds shoud have 3 columns. If they have more, there are manually added pararameters 
        if n_len_line_u > n_ids_and_functional:
            l_ids_and_functional         = line_unfilled[0:n_ids_and_functional]                                    #   1     2     1
            l_parameters_in_moleculetype = line_unfilled[n_ids_and_functional:n_len_line_u]                         # 0.1234       0.4321
            l_atom_names                 = line_filled[n_len_line_u:(n_len_line_u+n_ids_and_functional-1)]         # CA     CB

            ll_bonds_filled[i] = l_ids_and_functional + l_atom_names + l_parameters_in_moleculetype #reconstruct the line

    
    
    # [ angles ]
    for i in range(0, len(ll_angles)): #
        line_unfilled = ll_angles[i] #
        line_filled   = ll_angles_filled[i] #

        n_len_line_u = len(line_unfilled)
        n_len_line_f  = len(line_filled)

        #in [ angles ] Im looking for this type of filled columns in the current line
        #[ angles ]
        #;  ai    aj    ak funct            c0            c1            c2            c3
        #    2     1     3     5             0.1234       0.4321
        n_ids_and_functional = 4 # angles shoud have 4 columns. If they have more, there are manually added pararameters 
        if n_len_line_u > n_ids_and_functional:
            l_ids_and_functional         = line_unfilled[0:n_ids_and_functional]                                    #    2     1     3     5
            l_parameters_in_moleculetype = line_unfilled[n_ids_and_functional:n_len_line_u]                         # 0.1234       0.4321
            l_atom_names                 = line_filled[n_len_line_u:(n_len_line_u+n_ids_and_functional-1)]         # CA     CB      CD
            
            ll_angles_filled[i] = l_ids_and_functional + l_atom_names + l_parameters_in_moleculetype #reconstruct the line

     
    # [ dihedrals ] proper
    for i in range(0, len(ll_dihedrals_proper)): #
        line_unfilled = ll_dihedrals_proper[i] #
        line_filled   = ll_dihedrals_filled_proper[i] #

        n_len_line_u = len(line_unfilled)
        n_len_line_f  = len(line_filled)

        #in [ dihedrals ] Im looking for this type of filled columns in the current line
        #[ dihedrals ]
        #;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
        #    2     1     5     6     9            0.1234       0.4321
        n_ids_and_functional = 5 # dihedrals shoud have 5 columns. If they have more, there are manually added pararameters 
        if n_len_line_u > n_ids_and_functional:
            l_ids_and_functional         = line_unfilled[0:n_ids_and_functional]                                    #   2     1     5     6     9
            l_parameters_in_moleculetype = line_unfilled[n_ids_and_functional:n_len_line_u]                         # 0.1234       0.4321
            l_atom_names                 = line_filled[n_len_line_u:(n_len_line_u+n_ids_and_functional-1)]         # CA     CB      CD      CE
            
            ll_dihedrals_filled_proper[i] = l_ids_and_functional + l_atom_names + l_parameters_in_moleculetype #reconstruct the line

     
    # [ dihedrals ] improper
    for i in range(0, len(ll_dihedrals_improper)): #
        line_unfilled = ll_dihedrals_improper[i] #
        line_filled   = ll_dihedrals_filled_improper[i] #

        n_len_line_u = len(line_unfilled)
        n_len_line_f  = len(line_filled)

        #in [ dihedrals ] Im looking for this type of filled columns in the current line
        #[ dihedrals ]
        #;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
        #    2     1     5     6     9            0.1234       0.4321
        n_ids_and_functional = 5 # dihedrals shoud have 5 columns. If they have more, there are manually added pararameters 
        if n_len_line_u > n_ids_and_functional:
            l_ids_and_functional         = line_unfilled[0:n_ids_and_functional]                                    #   2     1     5     6     9
            l_parameters_in_moleculetype = line_unfilled[n_ids_and_functional:n_len_line_u]                         # 0.1234       0.4321
            l_atom_names                 = line_filled[n_len_line_u:(n_len_line_u+n_ids_and_functional-1)]         # CA     CB      CD      CE
            
            ll_dihedrals_filled_improper[i] = l_ids_and_functional + l_atom_names + l_parameters_in_moleculetype #reconstruct the line



    ############create list of tlc commands########## xxx
    print("CLEAN PIPE creating list of tcl commands")
    l_tcl_commads = [] #create list of lists to store all lines of a tcl script that will be sent in the end




    # generate a list of atom ids to write the local and global atom ids
    l_atom_ids = [row[0] for row in ll_atoms]
    #write local and global ids
    if l_atom_ids !=[]:
        print("CLEAN PIPE processing ids")


        #write local ids
        l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-local ids}}")
        l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
        l_tcl_commads.append("graphics top material Opaque")
        l_tcl_commads.append("graphics top color green")


        # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
        n_index_prev_mol = n_first_id -1
        for n_molecule_counter in range(1,n_molecules+1):

            for atom_id in l_atom_ids:

                coords_i = df_coordinates.loc[str(int(atom_id)+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()

                l_tcl_commads.append(f'graphics top text {{{coords_i.get("x")*10+0.3:.3f} {coords_i.get("y")*10+0.3:.3f} {coords_i.get("z")*10+0.3:.3f}}} "{str(atom_id)}" size 1')
                l_tcl_commads.append("display update")
                l_tcl_commads.append("mol off top")

            #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
            n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol

        #write global ids
        l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-global ids}}")
        l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
        l_tcl_commads.append("graphics top material Opaque")
        l_tcl_commads.append("graphics top color green")

        

        n_index_prev_mol = n_first_id -1
        for n_molecule_counter in range(1,n_molecules+1):

            for atom_id in l_atom_ids:

                coords_i = df_coordinates.loc[str(int(atom_id)+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()

                l_tcl_commads.append(f'graphics top text {{{coords_i.get("x")*10+0.3:.3f} {coords_i.get("y")*10+0.3:.3f} {coords_i.get("z")*10+0.3:.3f}}} "{str(int(atom_id)+n_index_prev_mol)}" size 1')
                l_tcl_commads.append("display update")
                l_tcl_commads.append("mol off top")


            #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
            n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol

    
    if "[ bonds ]" in what_to_plot and ll_bonds_filled !=[]:
        print("CLEAN PIPE processing [ bonds ]")
    
        #go throught the ids of columns that contain the functional, and the columns that contains the parameters
        count = 1 #this is to set the molecule name if its a parameter
        n_lenght = max(len(line) for line in ll_bonds_filled) #this is the lenght of the biggest line
        id_functional = 2
        id_fist_parameter = 5
        columns = [id_functional]+list(range(id_fist_parameter,n_lenght))#list containing the id of column with functional + the ids with parameters
        for column in columns: 
            
            #define name of molecule
            if column == id_functional:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-bonds-functional}}")
            else:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-bonds-c{count}}}")
                count = count + 1
    
            #define basic properties of the molecule
            l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
            l_tcl_commads.append("graphics top material Opaque")


            # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
            n_index_prev_mol = n_first_id -1
            for n_molecule_counter in range(1,n_molecules+1):
    
        
                for line in ll_bonds_filled:
                    
                    coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_j = df_coordinates.loc[str(int(line[1])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_m = algelin.calc_intermediate_point(coords_i,coords_j,50)
            
                    
                    try:
                        param = line[column] #get the column, that might not exist if there are less columns that than n_lenght
                        
                        #define color. but it will be white if the parameter is zero
                        if float(param) == 0:
                            l_tcl_commads.append("graphics top color white")
                        else:
                            l_tcl_commads.append("graphics top color blue")
                
                         #draw
                        l_tcl_commads.append(f"graphics top cylinder {{{coords_i.get('x')*10:.3f} {coords_i.get('y')*10:.3f} {coords_i.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} radius 0.1")
                        l_tcl_commads.append(f'graphics top text {{{coords_m.get("x")*10-0.3:.3f} {coords_m.get("y")*10-0.3:.3f} {coords_m.get("z")*10-0.3:.3f}}} "{param}" size 1')
                        l_tcl_commads.append("display update")
                        l_tcl_commads.append("mol off top")
                    except:
                        pass #if the line[column] is out of range 

                #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
                n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol
    
    
    if "[ angles ]" in what_to_plot and ll_angles_filled !=[]:
        print("CLEAN PIPE processing [ angles ]")
    
        #go throught the ids of columns that contain the functional, and the columns that contains the parameters
        count = 1 #this is to set the molecule name if its a parameter
        n_lenght = max(len(line) for line in ll_angles_filled) #this is the lenght of the biggest line
        id_functional = 3
        id_fist_parameter = 7
        columns = [id_functional]+list(range(id_fist_parameter,n_lenght))#list containing the id of column with functional + the ids with parameters
        for column in columns: 
            
            #define name of molecule
            if column == id_functional:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-angles-functional}}")
            else:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-angles-c{count}}}")
                count = count + 1
    
            #define basic properties of the molecule
            l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
            l_tcl_commads.append("graphics top material Opaque")
            


            # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
            n_index_prev_mol = n_first_id -1
            for n_molecule_counter in range(1,n_molecules+1):
            
                for line in ll_angles_filled:
                    
                    coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_j = df_coordinates.loc[str(int(line[1])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_k = df_coordinates.loc[str(int(line[2])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    
                    coords_i2 = algelin.calc_intermediate_point(coords_j,coords_i,25)
                    coords_k2 = algelin.calc_intermediate_point(coords_j,coords_k,25)
                    coords_m  = algelin.calc_intermediate_point(coords_i2,coords_k2,50)
            
                    
                    try:
                        param = line[column] #get the column, that might not exist if there are less columns that than n_lenght
                        
                        #define color. but it will be white if the parameter is zero
                        if float(param) == 0:
                            l_tcl_commads.append("graphics top color white")
                        else:
                            l_tcl_commads.append("graphics top color blue")
                
                         #draw
                        l_tcl_commads.append(f"graphics top triangle {{{coords_i2.get('x')*10:.3f} {coords_i2.get('y')*10:.3f} {coords_i2.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} {{{coords_k2.get('x')*10:.3f} {coords_k2.get('y')*10:.3f} {coords_k2.get('z')*10:.3f}}}")
                        l_tcl_commads.append(f'graphics top text {{{coords_m.get("x")*10:.3f} {coords_m.get("y")*10:.3f} {coords_m.get("z")*10:.3f}}} "{param}" size 1')
                        l_tcl_commads.append("display update")
                        l_tcl_commads.append("mol off top")
                    except:
                        pass #if the line[column] is out of range 
    
                #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
                n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol
    
    #proper dihedrals
    if "[ dihedrals ]" in what_to_plot and ll_dihedrals_filled_proper !=[]:
        print("CLEAN PIPE processing [ dihedrals ] all types but 2")
    
        #go throught the ids of columns that contain the functional, and the columns that contains the parameters
        count = 1 #this is to set the molecule name if its a parameter
        n_lenght = max(len(line) for line in ll_dihedrals_filled_proper) #this is the lenght of the biggest line
        id_functional = 4
        id_fist_parameter = 9
        columns = [id_functional]+list(range(id_fist_parameter,n_lenght))#list containing the id of column with functional + the ids with parameters
        for column in columns: 
            
            #define name of molecule
            if column == id_functional:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-dihedrals(prop)-functional}}")
            else:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-dihedrals(prop)-c{count}}}")
                count = count + 1
    
            #define basic properties of the molecule
            l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
            l_tcl_commads.append("graphics top material Opaque")



            # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
            n_index_prev_mol = n_first_id -1
            for n_molecule_counter in range(1,n_molecules+1):
        
                for line in ll_dihedrals_filled_proper:
                    
                    coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_j = df_coordinates.loc[str(int(line[1])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_k = df_coordinates.loc[str(int(line[2])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_l = df_coordinates.loc[str(int(line[3])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    
                    coords_m = algelin.calc_intermediate_point(coords_j,coords_k,50)
            
                    
                    try:
                        param = line[column] #get the column, that might not exist if there are less columns that than n_lenght
                        
                        #define color. but it will be white if the parameter is zero
                        if float(param) == 0:
                            l_tcl_commads.append("graphics top color white")
                        else:
                            l_tcl_commads.append("graphics top color blue")
                
                        #draw
                        l_tcl_commads.append(f"graphics top cone {{{coords_m.get('x')*10:.3f} {coords_m.get('y')*10:.3f} {coords_m.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} radius 0.2")
                        l_tcl_commads.append(f"graphics top cone {{{coords_m.get('x')*10:.3f} {coords_m.get('y')*10:.3f} {coords_m.get('z')*10:.3f}}} {{{coords_k.get('x')*10:.3f} {coords_k.get('y')*10:.3f} {coords_k.get('z')*10:.3f}}} radius 0.2")
                        l_tcl_commads.append(f'graphics top text {{{coords_m.get("x")*10+0.3:.3f} {coords_m.get("y")*10+0.3:.3f} {coords_m.get("z")*10+0.3:.3f}}} "{param}" size 1')
                        l_tcl_commads.append("display update")
                        l_tcl_commads.append("mol off top")
                    except:
                        pass #if the line[column] is out of range 
        
                #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
                n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol
        
    
    
    
    #improper dihedrals
    if "[ dihedrals ]" in what_to_plot and ll_dihedrals_filled_improper !=[]:
        print("CLEAN PIPE processing [ dihedrals ] type 2")
    
        #go throught the ids of columns that contain the functional, and the columns that contains the parameters
        count = 1 #this is to set the molecule name if its a parameter
        n_lenght = max(len(line) for line in ll_dihedrals_filled_improper) #this is the lenght of the biggest line
        id_functional = 4
        id_fist_parameter = 9
        columns = [id_functional]+list(range(id_fist_parameter,n_lenght))#list containing the id of column with functional + the ids with parameters
        for column in columns: 
            
            #define name of molecule
            if column == id_functional:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-dihedrals(impr)-functional}}")
            else:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-dihedrals(impr)-c{count}}}")
                count = count + 1
    
            #define basic properties of the molecule
            l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
            l_tcl_commads.append("graphics top material Transparent")
    

            # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
            n_index_prev_mol = n_first_id -1
            for n_molecule_counter in range(1,n_molecules+1):

            
                for line in ll_dihedrals_filled_improper:
                    
                    coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_j = df_coordinates.loc[str(int(line[1])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_k = df_coordinates.loc[str(int(line[2])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    coords_l = df_coordinates.loc[str(int(line[3])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    
                    try:
                        param = line[column] #get the column, that might not exist if there are less columns that than n_lenght
                        
                        #define color. but it will be white if the parameter is zero
                        if float(param) == 0:
                            l_tcl_commads.append("graphics top color white")
                        else:
                            l_tcl_commads.append("graphics top color green")
                
                        #draw
        
                        #lines around the external triangle
                        #l_tcl_commads.append(f"graphics top triangle {{{coords_l.get('x')*10:.3f} {coords_l.get('y')*10:.3f} {coords_l.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} {{{coords_k.get('x')*10:.3f} {coords_k.get('y')*10:.3f} {coords_k.get('z')*10:.3f}}}")
                        l_tcl_commads.append(f"graphics top line {{{coords_l.get('x')*10:.3f} {coords_l.get('y')*10:.3f} {coords_l.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}}")
                        l_tcl_commads.append(f"graphics top line {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} {{{coords_k.get('x')*10:.3f} {coords_k.get('y')*10:.3f} {coords_k.get('z')*10:.3f}}}")
                        l_tcl_commads.append(f"graphics top line {{{coords_k.get('x')*10:.3f} {coords_k.get('y')*10:.3f} {coords_k.get('z')*10:.3f}}} {{{coords_l.get('x')*10:.3f} {coords_l.get('y')*10:.3f} {coords_l.get('z')*10:.3f}}}")
        
                        #transparent triangle in the internal triangle
                        l_tcl_commads.append(f"graphics top triangle {{{coords_i.get('x')*10:.3f} {coords_i.get('y')*10:.3f} {coords_i.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} {{{coords_k.get('x')*10:.3f} {coords_k.get('y')*10:.3f} {coords_k.get('z')*10:.3f}}}")
                        l_tcl_commads.append(f'graphics top text {{{coords_i.get("x")*10+0.3:.3f} {coords_i.get("y")*10+0.3:.3f} {coords_i.get("z")*10+0.3:.3f}}} "{param}" size 1')
                        l_tcl_commads.append("display update")
                        l_tcl_commads.append("mol off top")
                    except:
                        pass #if the line[column] is out of range 
    
                #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
                n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol
    
    
    
    if "[ cmap ]" in what_to_plot and ll_cmap !=[]:
        print("CLEAN PIPE processing [ cmap ]")
    
        l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-cmap}}")
        l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
        l_tcl_commads.append("graphics top material Opaque")



        # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
        n_index_prev_mol = n_first_id -1
        for n_molecule_counter in range(1,n_molecules+1):
    
            for line in ll_cmap:
                
                coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                coords_j = df_coordinates.loc[str(int(line[1])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                coords_k = df_coordinates.loc[str(int(line[2])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                coords_l = df_coordinates.loc[str(int(line[3])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                coords_m = df_coordinates.loc[str(int(line[4])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
        
                coords_m1 = algelin.calc_intermediate_point(coords_j,coords_k,50)
                coords_m2 = algelin.calc_intermediate_point(coords_k,coords_l,50)
        
                coords_a = algelin.calc_intermediate_point(coords_j,coords_k,25)
                coords_b = algelin.calc_intermediate_point(coords_j,coords_k,75)
                coords_c = algelin.calc_intermediate_point(coords_k,coords_l,25)
                coords_d = algelin.calc_intermediate_point(coords_k,coords_l,75)
                
                l_tcl_commads.append("graphics top color purple")
                l_tcl_commads.append(f"graphics top cone {{{coords_m1.get('x')*10:.3f} {coords_m1.get('y')*10:.3f} {coords_m1.get('z')*10:.3f}}} {{{coords_a.get('x')*10:.3f} {coords_a.get('y')*10:.3f} {coords_a.get('z')*10:.3f}}} radius 0.3")
                l_tcl_commads.append(f"graphics top cone {{{coords_m1.get('x')*10:.3f} {coords_m1.get('y')*10:.3f} {coords_m1.get('z')*10:.3f}}} {{{coords_b.get('x')*10:.3f} {coords_b.get('y')*10:.3f} {coords_b.get('z')*10:.3f}}} radius 0.3")
                
                l_tcl_commads.append(f"graphics top cone {{{coords_m2.get('x')*10:.3f} {coords_m2.get('y')*10:.3f} {coords_m2.get('z')*10:.3f}}} {{{coords_c.get('x')*10:.3f} {coords_c.get('y')*10:.3f} {coords_c.get('z')*10:.3f}}} radius 0.3")
                l_tcl_commads.append(f"graphics top cone {{{coords_m2.get('x')*10:.3f} {coords_m2.get('y')*10:.3f} {coords_m2.get('z')*10:.3f}}} {{{coords_d.get('x')*10:.3f} {coords_d.get('y')*10:.3f} {coords_d.get('z')*10:.3f}}} radius 0.3")
                
                l_tcl_commads.append("display update")
                l_tcl_commads.append("mol off top")

            #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
            n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol
    
    
    if "[ pairs ]" in what_to_plot and ll_pairs !=[]:
        print("CLEAN PIPE processing [ pairs ]")
        
        l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-pairs}}")
        l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
        l_tcl_commads.append("graphics top color yellow")
        l_tcl_commads.append("graphics top material Opaque")


        # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
        n_index_prev_mol = n_first_id -1
        for n_molecule_counter in range(1,n_molecules+1):
        
            for line in ll_pairs:
                
                coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                coords_j = df_coordinates.loc[str(int(line[1])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                param = line[2]
                coords_m = algelin.calc_intermediate_point(coords_i,coords_j,50)
                
                l_tcl_commads.append(f"graphics top line {{{coords_i.get('x')*10:.3f} {coords_i.get('y')*10:.3f} {coords_i.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}}")
                l_tcl_commads.append(f'graphics top text {{{coords_m.get("x")*10-0.3:.3f} {coords_m.get("y")*10-0.3:.3f} {coords_m.get("z")*10-0.3:.3f}}} "{param}" size 1')
                l_tcl_commads.append("display update")
                l_tcl_commads.append("mol off top")

            #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
            n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol


    
    
    if "[ exclusions ]" in what_to_plot and ll_exclusions !=[]:
        print("CLEAN PIPE processing [ exclusions ]")
        
        l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-exclusion}}")
        l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
        l_tcl_commads.append("graphics top color orange")
        l_tcl_commads.append("graphics top material Transparent")


        # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
        n_index_prev_mol = n_first_id -1
        for n_molecule_counter in range(1,n_molecules+1):
        
            for line in ll_exclusions:
                
                #in exclusions, the first id ina line is a reference  atom, and all others in that line are exclusions with respect to that one
                coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                for second_atom_id in range(1,len(line)):
                    coords_j = df_coordinates.loc[str(int(line[second_atom_id])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
                    l_tcl_commads.append(f"graphics top cylinder {{{coords_i.get('x')*10:.3f} {coords_i.get('y')*10:.3f} {coords_i.get('z')*10:.3f}}} {{{coords_j.get('x')*10:.3f} {coords_j.get('y')*10:.3f} {coords_j.get('z')*10:.3f}}} radius 0.2")
                    l_tcl_commads.append("display update")
                    l_tcl_commads.append("mol off top")
    
            #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
            n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol
    
    
    if "[ position_restraints ]" in what_to_plot and ll_posres !=[]:
        print("CLEAN PIPE processing [ position_restraints ]")
        
        #go throught the ids of columns that contain the functional, and the columns that contains the parameters
        count = 1 #this is to set the molecule name if its a parameter
        n_lenght = max(len(line) for line in ll_posres) #this is the lenght of the biggest line
        id_functional = 1
        id_fist_parameter = 2
        columns = [id_functional]+list(range(id_fist_parameter,n_lenght))#list containing the id of column with functional + the ids with parameters
        for column in columns: 
            
            #define name of molecule
            if column == id_functional:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-posit...-functional}}")
            else:
                l_tcl_commads.append(f"mol load graphics {{{s_mol_name[0:5]}-posit...-c{count}}}")
                count = count + 1
    
            #define basic properties of the molecule
            l_tcl_commads.append("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
            l_tcl_commads.append("graphics top material Transparent")
            


            # go throught all the instantiations of molecules of a certain type that are present in the gro, appending plotting commands
            n_index_prev_mol = n_first_id -1
            for n_molecule_counter in range(1,n_molecules+1):
        
                for line in ll_posres:
                
                    coords_i = df_coordinates.loc[str(int(line[0])+n_index_prev_mol), ['x', 'y', 'z']].astype(float).to_dict()
    
                    try:
                        param = line[column] #get the column, that might not exist if there are less columns that than n_lenght
                        
                        #define color. but it will be white if the parameter is zero
                        if float(param) == 0:
                            l_tcl_commads.append("graphics top color white")
                        else:
                            l_tcl_commads.append("graphics top color red")
                
                        #draw
                        l_tcl_commads.append(f"graphics top sphere {{{coords_i.get('x')*10:.3f} {coords_i.get('y')*10:.3f} {coords_i.get('z')*10:.3f}}} radius 0.3")
                        l_tcl_commads.append(f'graphics top text {{{coords_i.get("x")*10+0.4:.3f} {coords_i.get("y")*10+0.4:.3f} {coords_i.get("z")*10+0.4:.3f}}} "{param}" size 1')
                        l_tcl_commads.append("display update")
                        l_tcl_commads.append("mol off top")
                    except:
                        pass #if the line[column] is out of range 


                #the for loop that goes trought all molecules of a certain type ends after this line that updates the n_index_prev_mol to be used in the next iteration
                n_index_prev_mol = n_index_prev_mol + n_atoms_in_mol



    #create temporary text file with all the lines of the tcl script
    with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as temp_file:
        # Join all strings with a newline and write them at once.
        temp_file.write("\n".join(l_tcl_commads))
    temp_file_path = temp_file.name
    print(f"CLEAN PIPE created temp file {temp_file_path} to store the tcl script that will be sent to vmd")
    
    send_command_to_vmd("source "+temp_file_path.replace("\\", "\\\\"))

    #clean temp files
    bricksFileSystem.delete(s_top_with_inclusions)
    bricksFileSystem.delete(temp_file_path)


def highlight_id(s_gro, n_id):
    """
    given an a global id, will draw a sphere around that atom in vmd
    the input id is te gromacs global id, as in the gro. I say that because the id numbering in vmd is different because it starts at 0 instead of 1. ignore the vmd numbering
    
    example:

    s_gro = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/3_NPT/npt.gro"
    cl.highlight_id(s_gro, 97)
    
    """

    #we get the coordinates of that id and draw an sphere there, so not to be misguided by the vmd id numbering, that is different
    coords = bricksGRO.coordinate_by_id(s_gro, n_id,'dictionary') #get the coordinates as a dictionary, for example {'x': 0.926, 'y': 1.383, 'z': 1.367}
    
    send_command_to_vmd(f"mol load graphics {{id {n_id} spherical highlight}}") #create a new molecule and set its name
    send_command_to_vmd("display resetview")#required after creating a new molecule so it doesnt have different coordinates and transformations
    send_command_to_vmd("graphics top material Transparent")
    send_command_to_vmd("graphics top color pink")
    send_command_to_vmd(f"graphics top sphere {{{coords.get('x')*10:.3f} {coords.get('y')*10:.3f} {coords.get('z')*10:.3f}}} radius 0.4")
    send_command_to_vmd("display update")





def see_forces(s_tpr_file,s_trr_file,force_threshold):

    """
    this function Il plot arrows in a trajectory file where the forcess are above the define threshhold

    this functions works by creating and sourcing a tcl script that updates the drawing of the errows every time the user mooves the frame slide on vmd


    example usage:
    
    s_tpr_file = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/4_PROD/prod.tpr"
    s_trr_file = r"//wsl$/Ubuntu/home/bioinformatician/MD/pepticat9_truss_in_water/4_PROD/prod.trr"
    force_threshold = 300  # in kJ/mol/nm

    cl.see_forces(s_tpr_file,s_trr_file,force_threshold)
    
    """


    # Load the GROMACS trajectory (.trr) and topology (.tpr)
    u = mda.Universe(s_tpr_file, s_trr_file, refresh_offsets=True)
    atoms = u.atoms

    # this a basic mda example on how to get information
    #u = mda.Universe(s_tpr_file, s_trr_file, refresh_offsets=True)
    #u.trajectory[0] # this set a certain frame directly in the object
    #a = u.atoms[123].position
    #b = u.atoms[123].force
    #c = u.atoms[123].name
    #d = u.atoms[123].resid

    # Prepare force data as a TCL dictionary
    force_data_tcl = "set force_data [dict create]\n"

    for frame_idx, ts in enumerate(u.trajectory):
        #in mdanalysus, when we chose a frame, like for example "u.trajectory[0]" or when we loop throught it (like it was done in this for loop),
        #the object is set to that frame, so that when we get corrdinates and forces of the atom, they will be from that specifc frame

    
        frame_data = []
        forces = atoms.forces 

        for i, force in enumerate(forces):

            force = np.array(force, dtype=float)  # Ensure force is a NumPy array.  example of a force np.array([150.0, -50.0, 200.0])


            # Compute the magnitude (length) of the force vector
            # Example calculation:
            # np.linalg.norm([150.0, -50.0, 200.0]) ≈ 250.0
            magnitude = np.linalg.norm(force) 
        
            if magnitude > force_threshold:
                #print(f"i {i}")
                #print(f"force {force}")
                #print(f"magnitude {magnitude}")
            
                pos = atoms[i].position 
                pos = np.array(pos, dtype=float)  # Ensure pos is a NumPy array.  example of a position: np.array([1.5, 2.0, 3.5])
                #print(f"pos {pos}")


                # Normalize the force vector (get its direction)
                # Example where force = [150, -50, 200] and magnitude = 250
                # direction = [150/250, -50/250, 200/250] ≈ [0.6, -0.2, 0.8]
                # If magnitude is very small (~0), set direction to zero to avoid division by zero
                direction = force / magnitude if magnitude > 1e-10 else np.zeros(3)

            
                # Compute the end position of the arrow (position + scaled direction)
                # Example:
                # pos = [1.5, 2.0, 3.5]
                # direction = [0.6, -0.2, 0.8]
                # Scaling by 2: direction * 2 = [1.2, -0.4, 1.6]
                # end_pos = [1.5, 2.0, 3.5] + [1.2, -0.4, 1.6] = [2.7, 1.6, 5.1]
                end_pos = pos + (direction * 2)


                # Append formatted force vector
                frame_data.append(
                    f"    {{{pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f} "
                    f"{end_pos[0]:.5f} {end_pos[1]:.5f} {end_pos[2]:.5f}}}"
                )

        # Only add frame data if it contains at least one force vector
        if frame_data:
            frame_entry = "\\\n".join(frame_data)  # Join forces on new lines with proper indentation
            force_data_tcl += f"dict set force_data {frame_idx} [list \\\n{frame_entry}]\n"

    # TCL script template
    tcl_script_template = f"""\
proc enable_force_visualization {{}} {{
    global vmd_frame force_data
    trace variable vmd_frame([molinfo top]) w draw_forces
}}

proc disable_force_visualization {{}} {{
    global vmd_frame
    trace vdelete vmd_frame([molinfo top]) w draw_forces
    draw delete all
}}

proc draw_arrow {{mol start end}} {{
    set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
    graphics $mol color red
    graphics $mol material Opaque
    graphics $mol cylinder $start $middle radius 0.2
    graphics $mol cone $middle $end radius 0.3
}}

proc draw_forces {{ name element op }} {{
    global vmd_frame force_data
    draw delete all
    set mol [molinfo top]

    if {{[dict exists $force_data $vmd_frame([molinfo top])]}} {{
        set forces [dict get $force_data $vmd_frame([molinfo top])]

        if {{[llength $forces] > 0}} {{
            foreach force_vector $forces {{
                set start [lrange $force_vector 0 2]
                set end [lrange $force_vector 3 5]
                draw_arrow $mol $start $end
            }}
        }}
    }}
}}

{force_data_tcl}




enable_force_visualization

"""



    #create temporary text file with the tcl script
    with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as temp_file:
        # Join all strings with a newline and write them at once.
        temp_file.write(tcl_script_template)
    temp_file_path = temp_file.name
    print(f"CLEAN PIPE created temp file {temp_file_path} to store the tcl script that will be sent to vmd")
    

    send_command_to_vmd("source "+temp_file_path.replace("\\", "\\\\"))



    #clean temp file
    #bricksFileSystem.delete(temp_file_path)

