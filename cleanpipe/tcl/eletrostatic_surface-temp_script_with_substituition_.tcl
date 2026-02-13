# If there is already a molecule, delete all its current representations
if {[molinfo num] > 0} {
    set mid [molinfo top]
    set num_reps [molinfo $mid get numreps]
    for {set i 0} {$i < $num_reps} {incr i} {
        mol delrep 0 $mid
    }
}

display projection orthographic

# Load density cube as a new molecule (this becomes top)
#mol new /data2/henrique/qm/cubes_hf/hf_Dt.cube type cube waitfor all
#mol new /data2/henrique/qm/dihedral_scan/ang000_homo39_Dt.cube type cube waitfor all
mol new /data2/henrique/qm/dihedral_scan/ang000_homo39_Dt.cube type cube waitfor all
set mid [molinfo top]

# Add ESP cube as an additional volumetric dataset to the SAME molecule
#mol addfile /data2/henrique/qm/cubes_hf/hf_ESP_clipped.cube type cube waitfor all
#mol addfile /data2/henrique/qm/dihedral_scan/ang000_homo39_ESP.cube type cube waitfor all
mol addfile /data2/henrique/qm/dihedral_scan/ang000_homo39_ESP.cube type cube waitfor all

display resetview

# --- Rep 0: Licorice (atoms) ---
mol representation Licorice 0.1 12.0 12.0
mol color ColorID 2
mol selection all
mol material Opaque
mol addrep $mid

# --- Rep 1: Isosurface from dataset 0, colored by dataset 1 ---
mol representation Isosurface 0.01 0 0 0 1 1
mol color Volume 1
mol selection all
mol material Transparent
mol addrep $mid











