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




# put a grey representation in all, that will be overiden my each specific type
mol representation Licorice 1.7 50 50
set text "all"
mol selection $text
mol material AOChalky
mol color ColorID 2 ;# grey
mol addrep top



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
