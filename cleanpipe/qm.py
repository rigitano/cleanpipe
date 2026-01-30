import psi4
import time
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio
import plotly.io as pio
from pathlib import Path
from PIL import Image
import os, math, glob
import plotly.graph_objects as go

from cleanpipe import bricksFileSystem



#THIS HERE WILL SAVE THE wfn_hf THAT CAME AS A RESULT OF psi4.optimize INTO CUBE FILES THAT i WILL BE ABLE TO VISUALIZE
"""
os.makedirs("cubes_hf_homo", exist_ok=True)

psi4.set_options({
    "cubeprop_tasks": ["DENSITY","ORBITALS", "ESP"],   # ESP gives ESP.cube (and Dt.cube), ORBITALS gives Psi_a_N.cube
    "cubeprop_orbitals": [homo],             # only the HOMO (alpha). For beta in UHF you'd use negative indices.

    "cubeprop_filepath": "cubes_hf_homo",
    "cubic_grid_spacing": [0.2, 0.2, 0.2],
    "cubic_grid_overage": [4.0, 4.0, 4.0],
})

psi4.cubeprop(wfn_hf)

# --- Rename for convenience ---
# Typical outputs you’ll see: ESP.cube, Dt.cube, Psi_a_<homo>.cube
for f in glob.glob("cubes_hf_homo/*.cube"):
    base = os.path.basename(f)
    new = os.path.join("cubes_hf_homo", f"hf_homo{homo}_{base}")
    os.rename(f, new)

print(f"Wrote HOMO cube for orbital index {homo} into cubes_hf_homo/")


"""





def qm(chosen_molecule, s_theory, s_basis, s_out_folder_name, b_optimize=True):
    """
    chosen_molecule:        molecule created using psi4.geometry
    s_theory:               a string containing "hf" or "dft" or "ccsd(t)"
    s_basis:                a string contining a basis set. ex: "6-31++G(d,p)"
    s_out_folder_name:      name of output to be created and where the cube files wil be saved (density, eletrostatic potential and wavefunction)
    b_optimize:             True if you want the position of the atoms to be optimized to avoid unrealistic molecules

    water = psi4.geometry('''
    O
    H 1 0.96
    H 1 0.96 2 104.5
    ''')
    energy, wfn = cl.qm(water, "hf", "6-31G(d)", "qm_hf_water", b_optimize=True)
    energy, wfn = cl.qm(water, "dft", "6-31++G(d,p)", "qm_dft_water", b_optimize=True)
    energy, wfn = cl.qm(water, "ccsd(t)", "cc-pVTZ", "qm_ccsd_water", b_optimize=True)

    """

    #its wise to clean eventual previous settings, otherwise settings of previous runs can be used unadvertly 
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()

    ######################################
    if s_theory == "hf":

        psi4.set_options({
            "reference": "rhf",
            "scf_type": "df",           # Density fitting for faster computations
            "e_convergence": 1e-6,      # Energy convergence threshold
            "d_convergence": 1e-6      # Density convergence threshold
        })

        # Perform the calculation of orbitals and energy
        if b_optimize: # optimize molecule geometry. the molecule is modified in place, the variable ame is just a pointer
            energy, wfn = psi4.optimize(f"hf/{s_basis}", molecule=chosen_molecule, return_wfn=True) 
        else: # molecule as it is, even if stretched or unrealistic
            energy, wfn = psi4.energy(f"hf/{s_basis}", molecule=chosen_molecule, return_wfn=True) 

        # Wait for 1 second to avoid bug in ps4, where the content of the variable with the energy has garbage on it if I try to print it too soon
        time.sleep(1)

    ############################################
    if s_theory == "dft":


        # Set calculation options
        psi4.set_options({
            "reference": "rks",
            "scf_type": "df",           # Density fitting for faster computations
            "e_convergence": 1e-6,      # Energy convergence threshold
            "d_convergence": 1e-6,      # Density convergence threshold
        })


        # Perform the calculation of orbitals and energy
        if b_optimize: # optimize molecule geometry. the molecule is modified in place, the variable ame is just a pointer
            energy, wfn = psi4.optimize(f"b3lyp/{s_basis}", molecule=chosen_molecule, return_wfn=True) 
        else: # molecule as it is, even if stretched or unrealistic
            energy, wfn = psi4.energy(f"b3lyp/{s_basis}", molecule=chosen_molecule, return_wfn=True) 



        # Wait for 1 second to avoid bug in ps4, where the content of the variable with the energy has garbage on it if I try to print it too soon
        time.sleep(1)

    #########################################
    if s_theory == "ccsd(t)":


        #ATENTION: you might want to do a geometry optimization using DFT first, to speed things up

        psi4.set_options({
            "reference": "rhf",   # Restricted Hartree-Fock
            "scf_type": "df",    # Use density fitting for SCF acceleration
            "cc_type": "df",
            
        })


        # Perform the calculation of orbitals and energy
        if b_optimize: # optimize molecule geometry. the molecule is modified in place, the variable ame is just a pointer
            energy, wfn = psi4.optimize(f"ccsd(t)/{s_basis}", molecule=chosen_molecule, return_wfn=True) 
        else: # molecule as it is, even if stretched or unrealistic
            energy, wfn = psi4.energy(f"ccsd(t)/{s_basis}", molecule=chosen_molecule, return_wfn=True) 



        # Wait for 1 second to avoid bug in ps4, where the content of the variable with the energy has garbage on it if I try to print it too soon
        time.sleep(1)

    #########################################
    else:
        print('the inserted theory is not supported. the possible choices are: "hf", "dft" or "ccsd(t)" ')
        return None
    #########################################
    #now lets export results, for hf and dft theories


    lumo_id = wfn.nalpha()
    homo_id = lumo_id - 1


    os.makedirs(s_out_folder_name, exist_ok=True)

    psi4.set_options({
        "cubeprop_tasks": ["DENSITY","ORBITALS", "ESP"],   # ESP gives ESP.cube (and Dt.cube), ORBITALS gives Psi_a_N.cube
        "cubeprop_orbitals": [homo_id,lumo_id],             # only the HOMO (alpha). For beta in UHF you'd use negative indices.

        "cubeprop_filepath": s_out_folder_name,
        "cubic_grid_spacing": [0.2, 0.2, 0.2],
        "cubic_grid_overage": [4.0, 4.0, 4.0],
    })

    psi4.cubeprop(wfn)

    # --- Rename for convenience ---
    # Typical outputs you’ll see: ESP.cube, Dt.cube, Psi_a_<homo>.cube
    #for f in glob.glob(f"{s_out_folder_name}/*.cube"):
    #    base = os.path.basename(f)
    #    new = os.path.join(s_out_folder_name, f"hf_homo{homo_id}_{base}")
    #    os.rename(f, new)

    print(f"CLEANPIPE MESSAGE cube files writen at {s_out_folder_name}")

    print(f"CLEANPIPE MESSAGE the energy is {energy}")

    print(energy)

    return energy, wfn

def qm_bond_scan(chosen_molecule, l_bond_indeces, s_theory='HF', s_basis='6-31G(d)', s_out_folder_name='bond_scan'):
    """

    cl.qm_bond_scan(cya, [0, 2], s_theory='HF', s_basis='6-31G(d)', s_out_folder_name='bond_scan')
    """

    bricksFileSystem.delete(s_out_folder_name)
    os.makedirs(s_out_folder_name, exist_ok=True)

    
    #results will be saved here if you need to see them after closing jupyter
    datafile = Path(s_out_folder_name) / "scan_data.txt"
   

    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()

    psi4.set_options({
        "reference": "rhf",
        "scf_type": "df",
        "e_convergence": 1e-9,
        "d_convergence": 1e-9,
    })

    mol0 = psi4.core.get_active_molecule().clone()
    energies = []

    for count, r in enumerate(np.arange(1.5, 5, 0.1)):
        mol = mol0.clone()

        geometric_keywords = {
            'coordsys': 'tric',
            'constraints': {
                'set': [{
                    'type': 'distance',
                    'indices': l_bond_indeces,
                    'value': float(r)
                }]
            }
        }

        try:
            E = psi4.optimize(
                f'{s_theory}/{s_basis}',
                molecule=mol,
                engine='geometric',
                optimizer_keywords=geometric_keywords
            )

            #store in variable
            energies.append((count, r, E))

            #store in data file, if you need the data after closing jupyter
            with open(datafile, "a") as f:
                f.write(f"{{'step': {count}, 'distance': {r}, 'energy': {E}}},\n")

            #save molecule geometry file
            mol.save_xyz_file(
                os.path.join(s_out_folder_name, f'optimized_distance_{r:.3f}.xyz'), 
                True
            )

            

        except Exception as e:
            with open(datafile, "a") as f:
                f.write(f"{{'step': {count}, 'distance': {r}, 'energy': ERROR}},\n")
            continue



def qm_angle_scan(chosen_molecule, l_angle_indeces, s_theory='HF', s_basis='6-31G(d)', s_out_folder_name='angle_scan'):
    """

    cl.qm_bond_scan(cya, [8, 0, 2], s_theory='HF', s_basis='6-31G(d)', s_out_folder_name='angle_scan')
    """

    bricksFileSystem.delete(s_out_folder_name)
    os.makedirs(s_out_folder_name, exist_ok=True)

    
    #results will be saved here if you need to see them after closing jupyter
    datafile = Path(s_out_folder_name) / "scan_data.txt"
   

    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()

    psi4.set_options({
        "reference": "rhf",
        "scf_type": "df",
        "e_convergence": 1e-9,
        "d_convergence": 1e-9,
    })

    mol0 = psi4.core.get_active_molecule().clone()
    energies = []

    for count, theta in enumerate(range(60, 181, 2)):
        mol = mol0.clone()

        geometric_keywords = {
            'coordsys': 'tric',
            'constraints': {
                'set': [{
                    'type': 'angle',
                    'indices': l_angle_indeces,
                    'value': float(theta)
                }]
            }
        }

        try:
            E = psi4.optimize(
                f'{s_theory}/{s_basis}',
                molecule=mol,
                engine='geometric',
                optimizer_keywords=geometric_keywords
            )

            #store in variable
            energies.append((count, theta, E))

            #store in data file, if you need the data after closing jupyter
            with open(datafile, "a") as f:
                f.write(f"{{'step': {count}, 'angle': {theta}, 'energy': {E}}},\n")

            #save molecule geometry file
            mol.save_xyz_file(
                os.path.join(s_out_folder_name, f'optimized_angle_{theta:03d}.xyz'), 
                True
            )

            

        except Exception as e:
            with open(datafile, "a") as f:
                f.write(f"{{'step': {count}, 'angle': {r}, 'energy': ERROR}},\n")
            continue





def qm_dihedral_scan(chosen_molecule, l_dihedral_indeces, s_theory='HF', s_basis='6-31G(d)', s_out_folder_name='dihedral_scan'):
    """

    cl.qm_dihedral_scan(cya, [8, 0, 2, 3], s_theory='HF', s_basis='6-31G(d)', s_out_folder_name='dihedral_scan')
    """

    bricksFileSystem.delete(s_out_folder_name)
    os.makedirs(s_out_folder_name, exist_ok=True)

    
    #results will be saved here if you need to see them after closing jupyter
    datafile = Path(s_out_folder_name) / "scan_data.txt"
   

    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()

    psi4.set_options({
        "reference": "rhf",
        "scf_type": "df",
        "e_convergence": 1e-9,
        "d_convergence": 1e-9,
    })

    mol0 = psi4.core.get_active_molecule().clone()
    energies = []

    for count, a in enumerate(range(0, 360, 5)): #id, final angle, step
        mol = mol0.clone()

        geometric_keywords = {
            'coordsys': 'tric',
            'constraints': {
                'set': [{
                    'type': 'dihedral',
                    'indices': l_dihedral_indeces,
                    'value': float(a)   # degrees
                }]
            }
        }

        try:
            E = psi4.optimize(
                f'{s_theory}/{s_basis}',
                molecule=mol,
                engine='geometric',
                optimizer_keywords=geometric_keywords
            )

            #store in variable
            energies.append((count, a, E))

            #store in data file, if you need the data after closing jupyter
            with open(datafile, "a") as f:
                f.write(f"{{'step': {count}, 'angle': {a}, 'energy': {E}}},\n")

            #save molecule geometry file
            mol.save_xyz_file(
                os.path.join(s_out_folder_name, f'optimized_torsion_dihedral_{a:03d}.xyz'), 
                True
            )

            

        except Exception as e:
            with open(datafile, "a") as f:
                f.write(f"{{'step': {count}, 'angle': {a}, 'energy': ERROR}},\n")
            continue





def read_cube(cube_file):
    """
    cube file reader to be used other functions
    """


    with open(cube_file, "r") as f:
        lines = f.readlines()

    natoms = int(lines[2].split()[0])
    origin = np.array([float(x) for x in lines[2].split()[1:4]], float)

    dims = []
    axes = []
    for i in range(3):
        p = lines[3+i].split()
        dims.append(int(p[0]))
        axes.append([float(p[1]), float(p[2]), float(p[3])])
    dims = tuple(dims)
    axes = np.array(axes, float)

    header_len = 2 + 1 + 3 + natoms
    toks = []
    for ln in lines[header_len:]:
        toks.extend(ln.split())
    data = np.array([float(x) for x in toks], float)

    vals = data.reshape(dims, order="C")
    spacing = np.array([np.linalg.norm(axes[0]), np.linalg.norm(axes[1]), np.linalg.norm(axes[2])], float)
    return origin, spacing, dims, vals



def read_cube_raw(cube_file):
    """
    another cube file reader to be used other functions
    """

    with open(cube_file, "r") as f:
        lines = f.readlines()

    natoms = int(lines[2].split()[0])
    nx = int(lines[3].split()[0])
    ny = int(lines[4].split()[0])
    nz = int(lines[5].split()[0])

    header_len = 2 + 1 + 3 + natoms
    toks = []
    for ln in lines[header_len:]:
        toks.extend(ln.split())
    data = np.array([float(x) for x in toks], dtype=float)
    return (nx, ny, nz), data




def cubes_check1(psi_cube_path):
    """
    you should chose the wavefunction cube file
    this will plot slices on a jupyter notebook so you can see this is a diagnistics to see if order should be "C" or "F".
    but this letter info should be properly set only in the reshape funcions

    example:
    cl.qm_check1("cubes_hf/Psi_a_38_38-A.cube"):
    """


    (dims, data) = read_cube_raw(psi_cube_path)
    nx, ny, nz = dims
    print("dims:", dims, "nvals:", data.size)

    vals_C = data.reshape((nx, ny, nz), order="C")
    vals_F = data.reshape((nx, ny, nz), order="F")

    ix, iy, iz = nx//2, ny//2, nz//2

    fig, axs = plt.subplots(2, 3, figsize=(12, 7))

    axs[0,0].imshow(vals_C[ix,:,:], origin="lower"); axs[0,0].set_title("C-order: x-mid")
    axs[0,1].imshow(vals_C[:,iy,:], origin="lower"); axs[0,1].set_title("C-order: y-mid")
    axs[0,2].imshow(vals_C[:,:,iz], origin="lower"); axs[0,2].set_title("C-order: z-mid")

    axs[1,0].imshow(vals_F[ix,:,:], origin="lower"); axs[1,0].set_title("F-order: x-mid")
    axs[1,1].imshow(vals_F[:,iy,:], origin="lower"); axs[1,1].set_title("F-order: y-mid")
    axs[1,2].imshow(vals_F[:,:,iz], origin="lower"); axs[1,2].set_title("F-order: z-mid")

    plt.tight_layout()
    plt.show()


def cubes_check2(density_cube,esp_cube):
    """
    you should chose the density cube and esp cube files
    this will plot the elctrostatic potential in a jupyter notebook

    example:
    cl.qm_check2("cubes_hf/Dt.cube", "cubes_hf/ESP.cube"):
    """

    origin, spacing, dims, rho = read_cube(density_cube)
    origin2, spacing2, dims2, esp = read_cube(esp_cube)

    # sanity check (must match)
    assert dims == dims2 and np.allclose(origin, origin2) and np.allclose(spacing, spacing2), \
        "Density and ESP cubes are not on the same grid."

    grid = pv.ImageData(dimensions=dims, spacing=spacing, origin=origin)

    # IMPORTANT: use Fortran order flatten for PyVista ImageData
    grid.point_data["rho"] = rho.ravel(order="F")
    grid.point_data["esp"] = esp.ravel(order="F")

    # classic molecular surface from density
    iso = 0.001
    surf = grid.contour([iso], scalars="rho")

    # Color that surface by ESP (sample/interpolate from grid onto the surface points)
    surf = surf.sample(grid)  # adds "esp" onto surf point data

    p = pv.Plotter(notebook=True)
    p.add_mesh(surf, scalars="esp", opacity=1.0)
    p.add_axes()
    p.show()

def cubes_check3(esp_cube):
    """
    yuo should chose the esp cube file
    this will print stats and slices on a jupyter notebook, so you can spot outliers in ESP that might mess up the automatic coloring in VMD
    if there are outliers, you should clip them out using the function clip_cube


    example:
    cl.cubes_check3("cubes_hf/ESP.cube"):
    """


    origin, spacing, dims, esp = read_cube(esp_cube)
    #origin, spacing, dims, rho = read_cube(density_cube)

    print("ESP dims:", dims, "spacing (bohr):", spacing)
    print("ESP min/max:", float(esp.min()), float(esp.max()))
    print("ESP mean/std:", float(esp.mean()), float(esp.std()))
    print("ESP fraction negative:", float((esp < 0).mean()))
    print("ESP fraction positive:", float((esp > 0).mean()))
    print("ESP fraction near zero (|v|<1e-3):", float((np.abs(esp) < 1e-3).mean()))

    # histogram (distribution)
    plt.figure(figsize=(6,4))
    plt.hist(esp.ravel(), bins=200)
    plt.xlabel("ESP value")
    plt.ylabel("count")
    plt.title("ESP cube histogram")
    plt.show()

    # plot 3 orthogonal mid-slices
    cx, cy, cz = [d//2 for d in dims]

    plt.figure(figsize=(6,5))
    plt.imshow(esp[:, :, cz].T, origin="lower")
    plt.colorbar(label="ESP")
    plt.title("ESP slice (z mid)")
    plt.show()

    plt.figure(figsize=(6,5))
    plt.imshow(esp[:, cy, :].T, origin="lower")
    plt.colorbar(label="ESP")
    plt.title("ESP slice (y mid)")
    plt.show()

    plt.figure(figsize=(6,5))
    plt.imshow(esp[cx, :, :].T, origin="lower")
    plt.colorbar(label="ESP")
    plt.title("ESP slice (x mid)")
    plt.show()



def clip_cube(infile, outfile, vmin=-0.2, vmax=0.2):
    """

    """

    with open(infile, "r") as f:
        lines = f.readlines()

    natoms = int(lines[2].split()[0])
    header_len = 2 + 1 + 3 + natoms

    header = lines[:header_len]
    data_tokens = []
    for ln in lines[header_len:]:
        data_tokens.extend(ln.split())

    data = np.array(data_tokens, dtype=float)
    data = np.clip(data, vmin, vmax)

    # Write back with 6 values per line (common cube formatting)
    with open(outfile, "w") as f:
        f.writelines(header)
        for i in range(0, data.size, 6):
            chunk = data[i:i+6]
            f.write(" ".join(f"{x: .6e}" for x in chunk) + "\n")




def plot_isosurfaces_from_expression(
    expr,
    *,
    L=20.0,
    n=64,
    what="both",          # "real", "imag", or "both"
    level=None,           # if None, auto from percentile of |component|
    percentile=99.0,      # used when level is None
    opacity=0.25,         # transparency so you can see through
    title=None,
    colors=("green", "yellow"),  # (real_color, imag_color)
):
    """
    Plot signed isosurfaces (positive/negative lobes) for Re(psi), Im(psi),
    or both overlaid.

    expr: string for psi(x,y,z). You may use x,y,z,r and numpy funcs like exp,sqrt,sin,...
          Use I for imaginary unit (1j).

    what:
      - "real": plot Re(psi) only
      - "imag": plot Im(psi) only
      - "both": overlay Re(psi) and Im(psi) (Re=green, Im=yellow by default)




    example:
    cl.plot_isosurfaces_from_expression("exp(-r)/sqrt(pi)",L=18, n=72)

    look at all the expressions that are solutions from the schrodinger equation for one electro hydrogen

    1s  exp(-r)/sqrt(pi)
    2s  (1/(4*sqrt(2*pi))) * (2 - r) * exp(-r/2)
    2pz (1/(4*sqrt(2*pi))) * z * exp(-r/2)
    2p  -(1/(8*sqrt(pi))) * (x + I*y) * exp(-r/2)
    2p  (1/(8*sqrt(pi))) * (x - I*y) * exp(-r/2)
    3s  (1/(81*sqrt(3*pi))) * (27 - 18*r + 2*r**2) * exp(-r/3)
    3   (1/(81*sqrt(6*pi))) * (2*z**2 - x**2 - y**2) * exp(-r/3)
    3   -(1/(81*sqrt(pi))) * z * (x + I*y) * exp(-r/3)
    3   (1/(81*sqrt(pi))) * z * (x - I*y) * exp(-r/3)
    3   (1/(162*sqrt(pi))) * (x + I*y)**2 * exp(-r/3)

    """
    

    # ---- grid ----
    x = np.linspace(-L, L, n)
    y = np.linspace(-L, L, n)
    z = np.linspace(-L, L, n)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)

    # ---- safe-ish eval namespace ----
    allowed = {name: getattr(np, name) for name in dir(np) if not name.startswith("_")}
    allowed.update({
        "x": X, "y": Y, "z": Z,
        "r": R, "R": R,
        "pi": np.pi,
        "I": 1j, "j": 1j
    })

    psi = eval(expr, {"__builtins__": {}}, allowed)

    fig = go.Figure()

    def add_component(component_name, V, color, level_override=None):
        # choose a symmetric level based on |V|
        absV = np.abs(V).ravel()
        lvl = level_override
        if lvl is None:
            lvl = np.percentile(absV, percentile)
        if lvl <= 0:
            # fallback, in case array is all zeros
            lvl = absV.max() if absV.max() > 0 else 1.0

        # single-color colorscale
        cs = [[0.0, color], [1.0, color]]

        # Positive lobe
        fig.add_trace(go.Isosurface(
            x=X.ravel(), y=Y.ravel(), z=Z.ravel(),
            value=V.ravel(),
            isomin=lvl, isomax=float(np.max(V)),
            surface_count=1,
            caps=dict(x_show=False, y_show=False, z_show=False),
            colorscale=cs,
            showscale=False,
            opacity=opacity,
            name=f"{component_name} +",
        ))

        # Negative lobe
        fig.add_trace(go.Isosurface(
            x=X.ravel(), y=Y.ravel(), z=Z.ravel(),
            value=V.ravel(),
            isomin=float(np.min(V)), isomax=-lvl,
            surface_count=1,
            caps=dict(x_show=False, y_show=False, z_show=False),
            colorscale=cs,
            showscale=False,
            opacity=opacity,
            name=f"{component_name} -",
        ))
        return lvl

    # ---- decide what to plot ----
    what = what.lower().strip()
    if what not in {"real", "imag", "both"}:
        raise ValueError("what must be 'real', 'imag', or 'both'")

    real_color, imag_color = colors
    used_levels = {}

    if what in {"real", "both"}:
        Vr = np.real(psi)
        used_levels["real"] = add_component("Re(ψ)", Vr, real_color, level_override=level)

    if what in {"imag", "both"}:
        Vi = np.imag(psi)
        # If user provided a single `level`, reuse it; otherwise auto separately for imag.
        level_imag = level if level is not None else None
        used_levels["imag"] = add_component("Im(ψ)", Vi, imag_color, level_override=level_imag)

    if title is None:
        if what == "both":
            title = f"Signed isosurfaces: Re(ψ) (green) + Im(ψ) (yellow)"
        else:
            title = f"Signed isosurfaces: {what.upper()}(ψ)"

    #this is a alternative that shows a grey grid behind with the axies
    """
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="x", yaxis_title="y", zaxis_title="z",
            aspectmode="cube",
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        legend=dict(itemsizing="constant"),
    )
    """

    fig.update_layout(
    title=title,
    scene=dict(
        aspectmode="cube",

        # Hide everything that looks like a box/grid
        xaxis=dict(visible=False, showbackground=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showbackground=False, showgrid=False, zeroline=False),
        zaxis=dict(visible=False, showbackground=False, showgrid=False, zeroline=False),
    ),
    margin=dict(l=0, r=0, t=40, b=0),
    legend=dict(itemsizing="constant"),
    paper_bgcolor="white",
    plot_bgcolor="white",
)
    fig.show()

    return used_levels  # handy to know the chosen thresholds





def phase_to_rgb_hsv(ph):
    """
    Map phase in [-pi, pi] to a cyclic RGB color wheel.
    Uses HSV with S=1, V=1 then converts to RGB (vectorized).
    """
    h = (ph + np.pi) / (2*np.pi)     # -> [0, 1]
    h = np.mod(h, 1.0)              # cyclic wrap

    i = np.floor(h * 6.0).astype(int)
    f = h * 6.0 - i
    i = i % 6

    p = np.zeros_like(h)
    q = 1.0 - f
    t = f

    r = np.choose(i, [1, q, p, p, t, 1])
    g = np.choose(i, [t, 1, 1, q, p, p])
    b = np.choose(i, [p, p, t, 1, 1, q])

    return np.stack([r, g, b], axis=1)

def read_cube_with_atoms(cube_file):
    with open(cube_file, "r") as f:
        lines = f.readlines()

    natoms = int(lines[2].split()[0])
    origin = np.array([float(x) for x in lines[2].split()[1:4]], float)  # bohr

    dims, axes = [], []
    for i in range(3):
        p = lines[3+i].split()
        dims.append(int(p[0]))
        axes.append([float(p[1]), float(p[2]), float(p[3])])  # bohr step vectors
    dims = tuple(dims)
    axes = np.array(axes, float)
    spacing = np.array([np.linalg.norm(axes[0]), np.linalg.norm(axes[1]), np.linalg.norm(axes[2])], float)

    # atom lines: atomic_number charge x y z  (coords in bohr)
    atom_lines = lines[6:6+natoms]
    atoms = []
    for ln in atom_lines:
        p = ln.split()
        Z = int(float(p[0]))
        x, y, z = map(float, p[2:5])
        atoms.append((Z, np.array([x, y, z], float)))
    atoms = [(Z, pos) for (Z, pos) in atoms]

    header_len = 2 + 1 + 3 + natoms
    toks = []
    for ln in lines[header_len:]:
        toks.extend(ln.split())
    data = np.array([float(x) for x in toks], float)

    psi = data.reshape(dims, order="C")  # your empirically correct reshape

    grid = pv.ImageData(dimensions=dims, spacing=spacing, origin=origin)
    grid.point_data["psi"] = psi.ravel(order="F")  # VTK point order

    return grid, psi, atoms

def build_molecule_mesh(atoms, bond_scale=1.2):
    """
    Create a PyVista mesh for atoms (spheres) and bonds (cylinders).
    Uses a simple distance-based bonding criterion using covalent radii.
    """
    # convert atom positions to Angstrom for radii convenience, then back to bohr for rendering consistency?
    # We'll render everything in bohr coordinates to match the cube grid.
    # So convert radii Å -> bohr when creating spheres/cylinders.
    ANG_TO_BOHR = 1.0 / BOHR_TO_ANG

    meshes = []

    # atom spheres
    for Z, pos_bohr in atoms:
        r_ang = COVALENT_RADIUS.get(Z, 0.7)
        r_bohr = r_ang * ANG_TO_BOHR
        sphere = pv.Sphere(radius=0.35*r_bohr, center=pos_bohr, theta_resolution=24, phi_resolution=24)
        meshes.append(("atom", Z, sphere))

    # bonds: naive O(N^2)
    n = len(atoms)
    for i in range(n):
        Zi, ri = atoms[i][0], COVALENT_RADIUS.get(atoms[i][0], 0.7)
        pi = atoms[i][1]
        for j in range(i+1, n):
            Zj, rj = atoms[j][0], COVALENT_RADIUS.get(atoms[j][0], 0.7)
            pj = atoms[j][1]
            d_bohr = np.linalg.norm(pj - pi)
            d_ang = d_bohr * BOHR_TO_ANG
            cutoff = bond_scale * (ri + rj)
            if d_ang < cutoff:
                # cylinder between pi and pj
                cyl_radius_bohr = 0.10 * ANG_TO_BOHR  # ~0.10 Å
                bond = pv.Cylinder(center=(pi+pj)/2, direction=(pj-pi), radius=cyl_radius_bohr, height=d_bohr, resolution=24)
                meshes.append(("bond", None, bond))

    return meshes




def create_frames_of_psi_waving(cube_path, out_dir, iso_value=0.05 ):
    """
    cube_path: you give the cube file that contains the wave function
    out_dir: folder name that wil be created to store the 180 frames of the video
    iso_value: you can also setup the isosurface value if you want the surface to be more far from the molecule

    this fucntion will create 180 video frames of lobes waving, always in oposit phases. all colors needed

    example:
    cl.create_frames_of_psi_waving("cubeprops_hf/Psi_a_38_38-A.cube", "psi_waving")
    """


    nframes=180


    GREEN  = np.array([0.0, 1.0, 0.0], float)
    YELLOW = np.array([1.0, 1.0, 0.0], float)

    # --- minimal element info ---
    # cube stores atomic number; we'll map a few common ones
    COVALENT_RADIUS = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66}  # Å (approx)
    ATOM_COLOR = {1: "white", 6: "gray", 7: "blue", 8: "red"}  # simple CPK-ish

    BOHR_TO_ANG = 0.529177210903


    # Load cube
    grid, psi, atoms = read_cube_with_atoms(cube_path)

    # fixed orbital magnitude isosurface
    grid.point_data["mag"] = np.abs(psi).ravel(order="F")
    surface = grid.contour([iso_value], scalars="mag")
    if surface.n_points == 0:
        raise RuntimeError("Isosurface empty. Try iso_value like 0.03.")

    # molecule geometry meshes
    mol_parts = build_molecule_mesh(atoms)

    # output folder
    os.makedirs(out_dir, exist_ok=True)
    for f in glob.glob(os.path.join(out_dir, "frame_*.png")):
        os.remove(f)

    # Render setup
    pv.global_theme.window_size = [1200, 912]
    plotter = pv.Plotter(off_screen=True)

    # add orbital surface (colored by rgb)
    surface.point_data["rgb"] = np.tile(GREEN, (surface.n_points, 1))
    plotter.add_mesh(surface, scalars="rgb", rgb=True, smooth_shading=True, opacity=0.75)

    # add molecule: atoms and bonds
    for kind, Z, mesh in mol_parts:
        if kind == "atom":
            color = ATOM_COLOR.get(Z, "lightgray")
            plotter.add_mesh(mesh, color=color, smooth_shading=True)
        else:
            plotter.add_mesh(mesh, color="darkgray", smooth_shading=True)

    plotter.camera_position = "iso"
    plotter.enable_anti_aliasing()
    plotter.set_background("white")

    plotter.reset_camera()

    cam = plotter.camera
    cam.Azimuth(50)      # rotate around vertical axis
    cam.Elevation(50)    # rotate up/down
    cam.Roll(40)         # tilt view

    plotter.render()     # apply the change once

    # Render frames
    for k in range(nframes):
        theta = 2.0 * math.pi * (k / nframes)
        c, s = math.cos(theta), math.sin(theta)

        re = psi * c
        im = -psi * s
        phase = np.arctan2(im, re)

        grid.point_data["phase"] = phase.ravel(order="F")
        sampled = surface.sample(grid)
        ph = np.array(sampled.point_data["phase"])

        rgb = phase_to_rgb_hsv(ph)
        surface.point_data["rgb"] = rgb


        plotter.render()
        plotter.screenshot(os.path.join(out_dir, f"frame_{k:04d}.png"))

    print(f"Wrote {nframes} frames to {out_dir}/")





    ####### after saving the frames. lets transform them into a gif video


    # Folder containing frames
    folder = Path(out_dir) 
    pattern = "frame_*.png"

    fps = 30
    frame_duration_ms = int(1000 / fps)

    files = sorted(folder.glob(pattern))
    assert files, f"No files found in {folder} matching {pattern}"

    frames = [Image.open(p).convert("RGB") for p in files]

    out_gif = folder / "psi_waving.gif"
    frames[0].save(
        out_gif,
        save_all=True,
        append_images=frames[1:],
        duration=frame_duration_ms,
        loop=0,
        optimize=True,
    )


#xxx not sure if this is ok
def see_labeled_molecule(cya):
    """
    the imput must be a molecule generated like so:

    water = psi4.geometry('''
    O
    H 1 0.96
    H 1 0.96 2 104.5
    ''')

    the output is a visualization in jupyter notebook of that molecule with labels over the atoms. 
    those labels will be suited just for the procedure of creating a set of dihedral angles

    """



    # Get XYZ string from Psi4 (different Psi4 versions expose different helpers)
    try:
        xyz = cya.to_string(dtype="xyz")
    except TypeError:
        xyz = cya.save_string_xyz()

    # ---- Visualization with atom IDs ----
    use_one_based_ids = True   # set False if you want 0-based IDs (Psi4-style)

    view = py3Dmol.view(width=700, height=520)
    view.addModel(xyz, "xyz")
    view.setStyle({"stick": {"radius": 0.18}, "sphere": {"scale": 0.28}})
    view.setBackgroundColor("0xFFFFFF")

    # Parse XYZ to place labels (works regardless of py3Dmol/3Dmol selection quirks)
    lines = [ln.strip() for ln in xyz.splitlines() if ln.strip()]
    nat = int(lines[0])
    atom_lines = lines[2:2+nat]  # skip nat + comment

    for i, ln in enumerate(atom_lines):
        sym, x, y, z = ln.split()[:4]
        idx = (i + 1) if use_one_based_ids else i
        label = f"{idx}:{sym}"

        view.addLabel(
            label,
            {
                "position": {"x": float(x), "y": float(y), "z": float(z)},
                "fontSize": 12,
                "fontColor": "black",
                "backgroundColor": "white",
                "backgroundOpacity": 0.6,
                "borderColor": "black",
                "borderOpacity": 0.3,
                "inFront": True
            }
        )

    view.zoomTo()
    view.show()







def make_string_gif(
    filename="string.gif",
    *,
    L=1.0,
    c=1.0,
    mode_n=1,
    A=1.0,
    phase=0.0,
    seconds=3.0,
    fps=30,
    nx=700,
    dpi=150,
    ylim=None,
    style="grey",   # "grey" or "bgy"
    line_width=3,
    bg="white"      # "white" or "transparent"
):
    """
    Generate a clean GIF of a string tied at both ends.

    style:
      - "grey": plain grey line
      - "bgy" : blue(+A) / grey(0) / yellow(-A) gradient along the string

    bg:
      - "white": white background
      - "transparent": transparent background (often nice for overlays)



    make_string_gif("string_grey3.gif", mode_n=3, A=1.0, seconds=2, fps=30, style="grey")
    make_string_gif("string_bgy3.gif",  mode_n=3, A=1.0, seconds=2, fps=30, style="bgy")
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize


    # --- Physics: fixed ends standing wave (normal mode) ---
    def fixed_ends_mode(x, t, n=1, A=1.0, L=1.0, c=1.0, phase=0.0):
        omega = n * np.pi * c / L
        return A * np.sin(n * np.pi * x / L) * np.cos(omega * t + phase)


    # --- Colormap: -A -> yellow, 0 -> grey, +A -> blue ---
    def yellow_grey_blue_cmap():
        colors = [
            (1.0, 1.0, 0.0),   # yellow
            (0.5, 0.5, 0.5),   # grey
            (0.0, 0.0, 1.0),   # blue
        ]
        return LinearSegmentedColormap.from_list("yellow_grey_blue", colors, N=256)


    def _make_linecollection(x, y, cmap, norm, lw=3):
        points = np.column_stack([x, y]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(y[:-1])          # color per segment from y
        lc.set_linewidth(lw)
        return lc




    x = np.linspace(0, L, nx)
    nframes = int(np.round(seconds * fps))
    t_values = np.linspace(0, seconds, nframes, endpoint=False)

    if ylim is None:
        ylim = 1.2 * abs(A)

    # Figure and axis: CLEAN (no axes, no ticks, no grid, no title)
    fig, ax = plt.subplots(figsize=(7, 3), dpi=dpi)
    ax.set_xlim(0, L)
    ax.set_ylim(-ylim, ylim)
    ax.axis("off")

    if bg == "transparent":
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    # Initial frame
    y0 = fixed_ends_mode(x, t_values[0], n=mode_n, A=A, L=L, c=c, phase=phase)

    artists = []

    if style.lower() == "grey":
        # Plain grey line (0.5 is mid-grey in Matplotlib)
        line, = ax.plot(x, y0, lw=line_width, color="0.5")
        artists = [line]

        def update(i):
            y = fixed_ends_mode(x, t_values[i], n=mode_n, A=A, L=L, c=c, phase=phase)
            line.set_ydata(y)
            return (line,)

    elif style.lower() in ("bgy", "blueyellow", "blue-yellow", "by"):
        cmap = yellow_grey_blue_cmap()
        norm = Normalize(vmin=-A, vmax=A)

        lc = _make_linecollection(x, y0, cmap=cmap, norm=norm, lw=line_width)
        ax.add_collection(lc)
        artists = [lc]

        def update(i):
            y = fixed_ends_mode(x, t_values[i], n=mode_n, A=A, L=L, c=c, phase=phase)

            # update geometry
            points = np.column_stack([x, y]).reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc.set_segments(segments)

            # update colors
            lc.set_array(y[:-1])
            return (lc,)

    else:
        plt.close(fig)
        raise ValueError('style must be "grey" or "bgy"')

    anim = FuncAnimation(fig, update, frames=nframes, blit=True)

    # Save GIF
    save_kwargs = {}
    if bg == "transparent":
        # Transparent background in GIF is supported with PillowWriter in many viewers
        save_kwargs["savefig_kwargs"] = {"transparent": True}

    anim.save(filename, writer=PillowWriter(fps=fps), dpi=dpi, **save_kwargs)
    plt.close(fig)
    return filename






def make_color_projected_gif(
    filename="color_projected.gif",
    *,
    L=1.0,
    c=1.0,
    mode_n=1,
    A=1.0,
    phase=0.0,
    seconds=3.0,
    fps=30,
    nx=700,
    dpi=150,
    style="bgy",     # "bgy" (recommended) or "grey"
    line_width=6,
    bg="white",      # "white" or "transparent"
    y0=0.0           # vertical position of the flat line
):
    """
    Creates a GIF where geometry is a flat horizontal line, but its color along x
    varies according to y(x,t) from a fixed-end standing wave.

    The "wave" is projected into color only: height does NOT change with time.

    make_color_projected_gif("projected_bgy3.gif", mode_n=3, A=1.0, style="bgy", seconds=2, fps=30)
    """


    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize


    def fixed_ends_mode(x, t, n=1, A=1.0, L=1.0, c=1.0, phase=0.0):
        omega = n * np.pi * c / L
        return A * np.sin(n * np.pi * x / L) * np.cos(omega * t + phase)


    def yellow_grey_blue_cmap():
        # -A -> yellow, 0 -> grey, +A -> blue
        colors = [
            (1.0, 1.0, 0.0),   # yellow
            (0.5, 0.5, 0.5),   # grey
            (0.0, 0.0, 1.0),   # blue
        ]
        return LinearSegmentedColormap.from_list("yellow_grey_blue", colors, N=256)



    x = np.linspace(0, L, nx)
    nframes = int(np.round(seconds * fps))
    t_values = np.linspace(0, seconds, nframes, endpoint=False)

    # Flat geometry: y is constant
    y_flat = np.full_like(x, float(y0))

    # Clean figure
    fig, ax = plt.subplots(figsize=(7, 1.2), dpi=dpi)
    ax.set_xlim(0, L)
    ax.set_ylim(y0 - 1.0, y0 + 1.0)  # just enough vertical room
    ax.axis("off")

    if bg == "transparent":
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    # Build segments for a flat line ONCE (geometry doesn't change)
    points = np.column_stack([x, y_flat]).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Initial color-signal from standing wave
    s0 = fixed_ends_mode(x, t_values[0], n=mode_n, A=A, L=L, c=c, phase=phase)

    if style.lower() == "grey":
        # Constant grey (no color undulation)
        lc = LineCollection(segments, colors=["0.5"], linewidths=line_width)
        ax.add_collection(lc)

        def update(i):
            return (lc,)

    elif style.lower() in ("bgy", "blueyellow", "blue-yellow", "by"):
        cmap = yellow_grey_blue_cmap()
        norm = Normalize(vmin=-A, vmax=A)

        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=line_width)
        lc.set_array(s0[:-1])  # color per segment from the signal
        ax.add_collection(lc)

        def update(i):
            s = fixed_ends_mode(x, t_values[i], n=mode_n, A=A, L=L, c=c, phase=phase)
            lc.set_array(s[:-1])  # update ONLY colors
            return (lc,)

    else:
        plt.close(fig)
        raise ValueError('style must be "grey" or "bgy"')

    anim = FuncAnimation(fig, update, frames=nframes, blit=True)

    save_kwargs = {}
    if bg == "transparent":
        save_kwargs["savefig_kwargs"] = {"transparent": True}

    anim.save(filename, writer=PillowWriter(fps=fps), dpi=dpi, **save_kwargs)
    plt.close(fig)
    return filename




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import hsv_to_rgb
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def fixed_ends_mode_complex(x, t, n=1, A=1.0, L=1.0, c=1.0, phase=0.0):
    """
    Complex standing-wave rotation:
        z(x,t) = A*sin(n*pi*x/L) * exp(i*(omega*t + phase))
    """
    omega = n * np.pi * c / L
    shape = A * np.sin(n * np.pi * x / L)
    return shape * np.exp(1j * (omega * t + phase))


def _phase_colors(z, sat=1.0, val=1.0):
    """
    Full color wheel from complex phase arg(z), with:
      +real (arg=0)  -> blue
      -real (arg=pi) -> yellow
    """
    arg = np.angle(z)  # [-pi, pi]
    hue = (2/3 + arg / (2 * np.pi)) % 1.0  # hue=2/3 is blue
    hsv = np.stack([hue, np.full_like(hue, sat), np.full_like(hue, val)], axis=-1)
    return hsv_to_rgb(hsv)

import numpy as np
from matplotlib.colors import hsv_to_rgb

def _phase_colors_with_amplitude_fade(z, A=1.0, sat_max=1.0, val=1.0, gamma=1.0):
    """
    Color by phase (full color wheel) but fade to grey near |z|=0 by reducing saturation.

    Requirements:
      +real -> blue, -real -> yellow (via hue shift)
      |z| ~ 0 -> grey-ish

    Parameters
    ----------
    z : complex array
    A : float
        Reference amplitude for normalization (same A as your wave).
    sat_max : float
        Saturation at full amplitude.
    val : float
        Brightness/value (keep 1.0 usually).
    gamma : float
        Controls how quickly color fades near zero:
          gamma > 1 fades more strongly near zero,
          gamma < 1 keeps more color near zero.
    """
    arg = np.angle(z)  # [-pi, pi]
    hue = (2/3 + arg / (2 * np.pi)) % 1.0  # +real -> blue

    # amplitude in [0,1]
    amp = np.abs(z) / max(A, 1e-12)
    amp = np.clip(amp, 0.0, 1.0)

    # Fade saturation to 0 near zero amplitude => grey
    sat = sat_max * (amp ** gamma)

    hsv = np.stack([hue, sat, np.full_like(hue, val)], axis=-1)
    return hsv_to_rgb(hsv)




def make_complex_string_gif_3d(
    filename="complex_string_3d.gif",
    *,
    L=1.0,
    c=1.0,
    mode_n=1,
    A=1.0,
    phase=0.0,
    seconds=3.0,
    fps=30,
    nx=700,
    dpi=150,
    style="phase",      # "grey" or "phase"
    line_width=3,
    bg="white",         # "white" or "transparent"
    elev=18,
    azim=-55,
    sat=1.0,
    val=1.0,
    axis_lw=1.0,        # thickness of the axis lines you asked for
):
    """
    3D view of complex standing wave rotating around the x-axis.

    Axes (SWITCHED as requested):
      - x: position along string
      - y: Im(z)  (horizontal complex axis)
      - z: Re(z)  (vertical real axis)  -> blue is up, yellow is down (via phase coloring)

    Visual style:
      - No grid, no ticks, no panes.
      - Draw only:
          * x-axis as a thin black line (y=z=0)
          * y-axis through origin (x=0, z=0)
          * z-axis through origin (x=0, y=0)


    make_complex_string_gif_3d("complex_axes_grey.gif",  mode_n=1, A=1.0, style="grey",  seconds=4, fps=30)
    make_complex_string_gif_3d("complex_axes_grey2.gif",  mode_n=2, A=1.0, style="grey",  seconds=4, fps=30)
    make_complex_string_gif_3d("complex_axes_grey3.gif",  mode_n=3, A=1.0, style="grey",  seconds=4, fps=30)
    make_complex_string_gif_3d("complex_axes_phase.gif", mode_n=1, A=1.0, style="phase", seconds=4, fps=30)
    make_complex_string_gif_3d("complex_axes_phase2.gif", mode_n=2, A=1.0, style="phase", seconds=4, fps=30)
    make_complex_string_gif_3d("complex_axes_phase3.gif", mode_n=3, A=1.0, style="phase", seconds=4, fps=30)
    """



    import numpy as np
    from matplotlib.colors import hsv_to_rgb

    def phase_color_with_grey_nodes(z, A=1.0, gamma=2.0, grey=0.5, sat=1.0, val=1.0):
        """
        Full color circle from phase, but enforce GREY (not white) near |z|=0.

        +Re (phase=0)   -> blue
        -Re (phase=pi)  -> yellow

        gamma controls how wide the grey node regions are.
        grey is the node color in RGB (0.5 = medium grey).
        """
        # hue from phase, shifted so +Re is blue
        arg = np.angle(z)
        hue = (2/3 + arg / (2 * np.pi)) % 1.0

        # vivid phase color
        hsv = np.stack([hue, np.full_like(hue, sat), np.full_like(hue, val)], axis=-1)
        rgb_phase = hsv_to_rgb(hsv)

        # amplitude weight: 0 near nodes -> 1 at antinodes
        amp = np.clip(np.abs(z) / max(A, 1e-12), 0.0, 1.0)
        w = amp ** gamma

        # blend grey <-> phase
        rgb_grey = np.full_like(rgb_phase, grey)
        rgb = (1 - w)[..., None] * rgb_grey + w[..., None] * rgb_phase
        return rgb


    def draw_box(ax, xlim, ylim, zlim, color="k", lw=1.0):
        x0, x1 = xlim
        y0, y1 = ylim
        z0, z1 = zlim

        corners = [
            (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
            (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
        ]

        edges = [
            (0,1),(1,2),(2,3),(3,0),  # bottom square
            (4,5),(5,6),(6,7),(7,4),  # top square
            (0,4),(1,5),(2,6),(3,7),  # vertical edges
        ]

        for i, j in edges:
            xi, yi, zi = corners[i]
            xj, yj, zj = corners[j]
            ax.plot([xi, xj], [yi, yj], [zi, zj], color=color, lw=lw)

    def draw_zy_square(ax, x_const, ylim, zlim, color="k", lw=1.0):
        """
        Draw ONLY the 4 edges of a square/rectangle in the zy-plane at x = x_const.

        Parameters
        ----------
        ax : 3D axes
        x_const : float
            The x position of the zy-plane (usually x0 = left boundary).
        ylim : (y0, y1)
            Limits along the y-axis (Im axis in your swapped setup).
        zlim : (z0, z1)
            Limits along the z-axis (Re axis in your swapped setup).
        """
        y0, y1 = ylim
        z0, z1 = zlim

        # 4 corners (x fixed)
        corners = [
            (x_const, y0, z0),
            (x_const, y1, z0),
            (x_const, y1, z1),
            (x_const, y0, z1),
        ]

        # connect edges: 0-1-2-3-0
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

        for i, j in edges:
            xi, yi, zi = corners[i]
            xj, yj, zj = corners[j]
            ax.plot([xi, xj], [yi, yj], [zi, zj], color=color, lw=lw)

    
    x = np.linspace(0, L, nx)
    nframes = int(np.round(seconds * fps))
    t_values = np.linspace(0, seconds, nframes, endpoint=False)

    # Initial curve
    zz0 = fixed_ends_mode_complex(x, t_values[0], n=mode_n, A=A, L=L, c=c, phase=phase)
    y0 = np.imag(zz0)   # complex axis (horizontal)
    z0 = np.real(zz0)   # real axis (vertical)

    fig = plt.figure(figsize=(7, 4), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")

    # Background
    if bg == "transparent":
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    # Limits (fixed)
    lim = 1.2 * abs(A)
    ax.set_xlim(0, L)
    ax.set_ylim(-lim, lim)  # Im
    ax.set_zlim(-lim, lim)  # Re (vertical)

    draw_zy_square(ax, x_const=0, ylim=(-lim, lim), zlim=(-lim, lim), color="k", lw=1)


    # Camera
    ax.view_init(elev=elev, azim=azim)

    # Remove all default axis clutter (grid, ticks, panes, etc.)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    # Hide axis panes completely
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            pass
            #axis.pane.set_facecolor((1, 1, 1, 0))
            #axis.pane.set_edgecolor((1, 1, 1, 0))
        except Exception:
            pass

    # Also hide the default axis lines (we'll draw our own)
    try:
        ax.w_xaxis.line.set_color((1, 1, 1, 0))
        ax.w_yaxis.line.set_color((1, 1, 1, 0))
        ax.w_zaxis.line.set_color((1, 1, 1, 0))
    except Exception:
        pass

    # --- Draw ONLY the three axes crossing the origin ---
    # x-axis: y=z=0, from 0..L (thin black line)
    ax.plot([0, L], [0, 0], [0, 0], color="k", lw=axis_lw)

    # y-axis (Im): x=0, z=0, from -lim..lim
    ax.plot([0, 0], [-lim, lim], [0, 0], color="k", lw=axis_lw)

    # z-axis (Re): x=0, y=0, from -lim..lim
    ax.plot([0, 0], [0, 0], [-lim, lim], color="k", lw=axis_lw)

    # --- Animated wave as a 3D segmented line (for per-segment colors) ---
    pts0 = np.column_stack([x, y0, z0])
    segs0 = np.stack([pts0[:-1], pts0[1:]], axis=1)
    lc = Line3DCollection(segs0, linewidths=line_width)
    ax.add_collection3d(lc)

    if style.lower() == "grey":
        lc.set_color((0.5, 0.5, 0.5, 1.0))

        def update(i):
            zz = fixed_ends_mode_complex(x, t_values[i], n=mode_n, A=A, L=L, c=c, phase=phase)
            yy = np.imag(zz)
            zz_re = np.real(zz)

            pts = np.column_stack([x, yy, zz_re])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            lc.set_segments(segs)
            return (lc,)

    elif style.lower() in ("phase", "colored", "color", "hsv"):
        #rgb = _phase_colors(zz0, sat=sat, val=val)
        rgb = phase_color_with_grey_nodes(zz0, A=A, gamma=2.0, grey=0.5, sat=sat, val=val)

        lc.set_color(np.column_stack([rgb[:-1], np.ones(nx - 1)]))

        def update(i):
            zz = fixed_ends_mode_complex(x, t_values[i], n=mode_n, A=A, L=L, c=c, phase=phase)
            yy = np.imag(zz)
            zz_re = np.real(zz)

            pts = np.column_stack([x, yy, zz_re])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            lc.set_segments(segs)

            #rgb = _phase_colors(zz, sat=sat, val=val)
            rgb = _phase_colors_with_amplitude_fade(zz, A=A, sat_max=sat, val=val, gamma=2.0)

            lc.set_color(np.column_stack([rgb[:-1], np.ones(nx - 1)]))
            return (lc,)

    else:
        plt.close(fig)
        raise ValueError('style must be "grey" or "phase"')

    anim = FuncAnimation(fig, update, frames=nframes, blit=True)

    save_kwargs = {}
    if bg == "transparent":
        save_kwargs["savefig_kwargs"] = {"transparent": True}




    anim.save(filename, writer=PillowWriter(fps=fps), dpi=dpi, **save_kwargs)
    plt.close(fig)
    return filename


def make_complex_color_rod_gif_3d(
    filename="complex_color_rod_3d.gif",
    *,
    L=1.0,
    c=1.0,
    mode_n=1,
    A=1.0,
    phase=0.0,
    seconds=3.0,
    fps=30,
    nx=700,
    dpi=150,
    style="phase",      # "grey" or "phase"
    rod_width=8,        # thickness of the rod
    bg="white",         # "white" or "transparent"
    elev=18,
    azim=-55,
    sat=1.0,
    val=1.0,
    axis_lw=1.0,
    zy_frame=True,      # draw the zy square at x=0
    x_axis_line=True,   # draw thin x-axis line under the rod
):
    """
    3D "rod" along the x-axis (y=z=0) whose color varies along x in time,
    using the SAME complex wave logic as the rotating curve.

    Axes convention (same as before, swapped):
      x: position
      y: Im(z)  (horizontal complex axis)
      z: Re(z)  (vertical real axis)

    But the geometry shown is ONLY the rod on the x-axis. The wave influences color only.


    make_complex_color_rod_gif_3d("rod_phase.gif", mode_n=1, A=1.0, style="phase", seconds=4, fps=30, rod_width=10)
    make_complex_color_rod_gif_3d("rod_phase2.gif", mode_n=2, A=1.0, style="phase", seconds=4, fps=30, rod_width=10)
    make_complex_color_rod_gif_3d("rod_phase3.gif", mode_n=3, A=1.0, style="phase", seconds=4, fps=30, rod_width=10)
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.colors import hsv_to_rgb
    from mpl_toolkits.mplot3d.art3d import Line3DCollection


    def fixed_ends_mode_complex(x, t, n=1, A=1.0, L=1.0, c=1.0, phase=0.0):
        omega = n * np.pi * c / L
        shape = A * np.sin(n * np.pi * x / L)
        return shape * np.exp(1j * (omega * t + phase))


    def _phase_colors(z, sat=1.0, val=1.0):
        """
        Phase -> full color wheel, with:
        +real (arg=0)  -> blue
        -real (arg=pi) -> yellow
        """
        arg = np.angle(z)
        hue = (2/3 + arg / (2 * np.pi)) % 1.0
        hsv = np.stack([hue, np.full_like(hue, sat), np.full_like(hue, val)], axis=-1)
        return hsv_to_rgb(hsv)

    import numpy as np
    from matplotlib.colors import hsv_to_rgb

    def _phase_colors_with_amplitude_fade(z, A=1.0, sat_max=1.0, val=1.0, gamma=1.0):
        """
        Color by phase (full color wheel) but fade to grey near |z|=0 by reducing saturation.

        Requirements:
        +real -> blue, -real -> yellow (via hue shift)
        |z| ~ 0 -> grey-ish

        Parameters
        ----------
        z : complex array
        A : float
            Reference amplitude for normalization (same A as your wave).
        sat_max : float
            Saturation at full amplitude.
        val : float
            Brightness/value (keep 1.0 usually).
        gamma : float
            Controls how quickly color fades near zero:
            gamma > 1 fades more strongly near zero,
            gamma < 1 keeps more color near zero.
        """
        arg = np.angle(z)  # [-pi, pi]
        hue = (2/3 + arg / (2 * np.pi)) % 1.0  # +real -> blue

        # amplitude in [0,1]
        amp = np.abs(z) / max(A, 1e-12)
        amp = np.clip(amp, 0.0, 1.0)

        # Fade saturation to 0 near zero amplitude => grey
        sat = sat_max * (amp ** gamma)

        hsv = np.stack([hue, sat, np.full_like(hue, val)], axis=-1)
        return hsv_to_rgb(hsv)


    def draw_zy_square(ax, x_const, ylim, zlim, color="k", lw=1.0):
        y0, y1 = ylim
        z0, z1 = zlim
        corners = [
            (x_const, y0, z0),
            (x_const, y1, z0),
            (x_const, y1, z1),
            (x_const, y0, z1),
        ]
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        for i, j in edges:
            xi, yi, zi = corners[i]
            xj, yj, zj = corners[j]
            ax.plot([xi, xj], [yi, yj], [zi, zj], color=color, lw=lw)


    import numpy as np
    from matplotlib.colors import hsv_to_rgb

    def phase_color_with_grey_nodes(z, A=1.0, gamma=2.0, grey=0.5, sat=1.0, val=1.0):
        """
        Full color circle from phase, but enforce GREY (not white) near |z|=0.

        +Re (phase=0)   -> blue
        -Re (phase=pi)  -> yellow

        gamma controls how wide the grey node regions are.
        grey is the node color in RGB (0.5 = medium grey).
        """
        # hue from phase, shifted so +Re is blue
        arg = np.angle(z)
        hue = (2/3 + arg / (2 * np.pi)) % 1.0

        # vivid phase color
        hsv = np.stack([hue, np.full_like(hue, sat), np.full_like(hue, val)], axis=-1)
        rgb_phase = hsv_to_rgb(hsv)

        # amplitude weight: 0 near nodes -> 1 at antinodes
        amp = np.clip(np.abs(z) / max(A, 1e-12), 0.0, 1.0)
        w = amp ** gamma

        # blend grey <-> phase
        rgb_grey = np.full_like(rgb_phase, grey)
        rgb = (1 - w)[..., None] * rgb_grey + w[..., None] * rgb_phase
        return rgb


    x = np.linspace(0, L, nx)
    nframes = int(np.round(seconds * fps))
    t_values = np.linspace(0, seconds, nframes, endpoint=False)

    lim = 1.2 * abs(A)

    fig = plt.figure(figsize=(7, 4), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")

    if bg == "transparent":
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    ax.set_xlim(0, L)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    ax.view_init(elev=elev, azim=azim)

    # Remove all default clutter
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            axis.pane.set_facecolor((1, 1, 1, 0))
            axis.pane.set_edgecolor((1, 1, 1, 0))
        except Exception:
            pass
    try:
        ax.w_xaxis.line.set_color((1, 1, 1, 0))
        ax.w_yaxis.line.set_color((1, 1, 1, 0))
        ax.w_zaxis.line.set_color((1, 1, 1, 0))
    except Exception:
        pass

    # Optional: draw the minimal reference axes you wanted
    if x_axis_line:
        ax.plot([0, L], [0, 0], [0, 0], color="k", lw=axis_lw)
    ax.plot([0, 0], [-lim, lim], [0, 0], color="k", lw=axis_lw)   # y-axis through origin
    ax.plot([0, 0], [0, 0], [-lim, lim], color="k", lw=axis_lw)   # z-axis through origin

    if zy_frame:
        draw_zy_square(ax, x_const=0, ylim=(-lim, lim), zlim=(-lim, lim), color="k", lw=axis_lw)

    # --- Build the rod segments ONCE (geometry fixed on x-axis) ---
    y_rod = np.zeros_like(x)
    z_rod = np.zeros_like(x)
    pts = np.column_stack([x, y_rod, z_rod])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)

    rod = Line3DCollection(segs, linewidths=rod_width)
    ax.add_collection3d(rod)

    # Initial color signal from complex wave
    zc0 = fixed_ends_mode_complex(x, t_values[0], n=mode_n, A=A, L=L, c=c, phase=phase)

    if style.lower() == "grey":
        rod.set_color((0.5, 0.5, 0.5, 1.0))

        def update(i):
            return (rod,)

    elif style.lower() in ("phase", "colored", "color", "hsv"):
        #rgb0 = _phase_colors(zc0, sat=sat, val=val)
        rgb0 = phase_color_with_grey_nodes(zc0, A=A, gamma=2.0, grey=0.5, sat=sat, val=val)

        rod.set_color(np.column_stack([rgb0[:-1], np.ones(nx - 1)]))

        def update(i):
            zc = fixed_ends_mode_complex(x, t_values[i], n=mode_n, A=A, L=L, c=c, phase=phase)
            rgb = rgb = _phase_colors_with_amplitude_fade(zc, A=A, sat_max=sat, val=val, gamma=2.0)
            rod.set_color(np.column_stack([rgb[:-1], np.ones(nx - 1)]))
            return (rod,)

    else:
        plt.close(fig)
        raise ValueError('style must be "grey" or "phase"')

    anim = FuncAnimation(fig, update, frames=nframes, blit=True)

    save_kwargs = {}
    if bg == "transparent":
        save_kwargs["savefig_kwargs"] = {"transparent": True}

    anim.save(filename, writer=PillowWriter(fps=fps), dpi=dpi, **save_kwargs)
    plt.close(fig)
    return filename

