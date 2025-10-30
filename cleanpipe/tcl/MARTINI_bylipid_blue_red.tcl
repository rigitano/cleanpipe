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





# reddish tone
mol representation Licorice 2.0 12.0 12.0
mol color ColorID 29
mol selection {resname SED_WILL_REPLACE_THIS}
mol material AOChalky
mol addrep top

# bluish tone
mol representation Licorice 2.0 12.0 12.0
mol color ColorID 23
mol selection {not resname SED_WILL_REPLACE_THIS}
mol material AOChalky
mol addrep top






# PROTEIN

set text "name BB"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 2 ;# grey
mol addrep top

mol selection {name regexp "^SC.*"}
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 2 ;# grey
mol addrep top