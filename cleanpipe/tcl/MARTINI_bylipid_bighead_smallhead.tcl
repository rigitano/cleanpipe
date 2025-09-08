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





mol representation Licorice 2.0 12.0 12.0
mol color ColorID 25
mol selection {resname regexp ".*(PC|SM|PG)$"}
mol material AOChalky
mol addrep top


mol representation Licorice 2.0 12.0 12.0
mol color ColorID 13
mol selection {resname regexp ".*(PE|PA|PS|CL)$"}
mol material AOChalky
mol addrep top

