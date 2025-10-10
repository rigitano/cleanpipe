#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}
    

# See protein as licorice
mol selection {protein}
mol representation Licorice 0.3 12 12
mol color Type ;# by atom type
mol material Goodsell
mol addrep top

# Add a cartoon representation for the backbone
mol selection {protein}
mol representation NewCartoon
mol color ColorID 2 ;# grey
mol material Goodsell
mol addrep top

# Visualize molecules within hydrogen bond range of the protein

mol selection {not protein and within 3.5 of protein}
mol representation Licorice 0.3 12 12
mol color Name
mol material Opaque
mol addrep top

# Visualize hydrogen bonds using Hbonds drawing method
mol selection {protein or (not protein and within 3.5 of protein)}
mol representation Hbonds
mol color ColorID 27 ;# magenta
mol material Opaque
mol addrep top