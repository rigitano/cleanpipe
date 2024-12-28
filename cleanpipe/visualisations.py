import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from matplotlib.ticker import ScalarFormatter, LogLocator
import seaborn as sns
import pandas as pd
import numpy as np



def plot_hbonds(s_ss_file,s_ps_file,s_pp_file,s_subtitle=""):
    """
    the inputs are 3 xvg files that have to be generated calling the tool 'gmx hbonds' 3 times:
    one to calculate solvent-solvent hbonds, other to calculate protein-solvent hbonds, and other to calculate protein-protein hbons

    the last option input is a subtitle

    example:
    cl.plot_hbonds("hb_ss.xvg","hb_ps.xvg","hb_pp.xvg","hello")
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

plot_ramachandran("rama.csv","hello")

