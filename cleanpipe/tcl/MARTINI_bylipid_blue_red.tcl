#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}



display shadows on
display ambientocclusion on
display aoambient 1.0
display aodirect 0.1
display rendermode GLSL




# Beautiful DBPC (bluish) and DLPC (reddish) visualization

# DBPC - bluish tone
mol representation Licorice 2.0 12.0 12.0
mol color ColorID 23
mol selection {resname DBPC}
mol material AOChalky
mol addrep top

# DLPC - reddish tone
mol representation Licorice 2.0 12.0 12.0
mol color ColorID 29
mol selection {resname DLPC}
mol material AOChalky
mol addrep top
