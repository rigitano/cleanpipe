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