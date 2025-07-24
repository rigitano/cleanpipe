#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}
    

# See protein as licorice
mol selection {protein}
mol representation Licorice
mol color Type ;# by atom type
mol material Opaque
mol addrep top

# Add a cartoon representation for the backbone
mol selection {protein}
mol representation cartoon
mol color ColorID 2 ;# grey
mol material Opaque
mol addrep top