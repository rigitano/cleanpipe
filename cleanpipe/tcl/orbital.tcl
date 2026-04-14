# If there is already a molecule, delete all its current representations
if {[molinfo num] > 0} {
    set mid [molinfo top]
    set num_reps [molinfo $mid get numreps]
    for {set i 0} {$i < $num_reps} {incr i} {
        mol delrep 0 $mid
    }
}

display projection orthographic




# ---- ORBITAL (HOMO/LUMO) ISOSURFACES: + blue, - yellow ----
# Add orbital cube as another volumetric dataset on the SAME molecule
#mol addfile /data2/henrique/qm/cubes_hf_homo/hf_homo38_hf_homo38_Psi_a_38_38-A.cube type cube waitfor all
#mol addfile /data2/henrique/qm/dihedral_scan/ang000_homo39_Psi_a_39_39-A.cube type cube waitfor all
#mol new SED_WILL_REPLACE_THIS_1 type cube waitfor all
set mid [molinfo top]

display resetview

# --- Rep 0: Licorice (atoms) ---
mol representation Licorice 0.1 12.0 12.0
mol color ColorID 2
mol selection all
mol material Opaque
mol addrep $mid


set orb_vol [expr {[molinfo $mid get numvolumedata] - 1}]

# Positive lobe (+iso): blue
mol representation Isosurface 0.01 $orb_vol 0 0 1 1
mol color ColorID 0
mol selection all
mol material Transparent
mol addrep $mid
set rep_pos [expr {[molinfo $mid get numreps] - 1}]


# Negative lobe (-iso): yellow
mol representation Isosurface -0.01 $orb_vol 0 0 1 1
mol color ColorID 4
mol selection all
mol material Transparent
mol addrep $mid
set rep_neg [expr {[molinfo $mid get numreps] - 1}]






