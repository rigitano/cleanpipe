#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}

display projection orthographic  

# See protein as licorice
mol selection {not water}
mol representation Licorice 0.2 12 12
mol color Type ;# by atom type
mol material Opaque ;# before, it was Goodsell.
mol addrep top

# Add a cartoon representation for the backbone
mol selection {protein}
mol representation NewCartoon
mol color ColorID 2 ;# grey
mol material Opaque ;# before, it was Goodsell
mol addrep top

# Visualize molecules within hydrogen bond range of the protein

mol selection {water and within 3.5 of (not water)}
mol representation Licorice 0.2 12 12
mol color Name
mol material Opaque
mol addrep top

# Visualize hydrogen bonds using Hbonds drawing method (between protein and water)
mol selection {(not water) or (water and within 3.5 of (not water))}
mol representation Hbonds
mol color ColorID 27 ;# magenta
mol material Opaque
mol addrep top

# Visualize hydrogen bonds using Hbonds drawing method (between atoms on the protein backbone)
mol selection {name N H CO O}
mol representation Hbonds
mol color ColorID 27 ;# magenta
mol material Opaque
mol addrep top