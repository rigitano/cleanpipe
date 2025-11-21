# Delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}


display projection orthographic

# first Ill set everithing as grey

set text "not name W"
mol selection $text
mol representation Licorice 0.3 12 12
mol material AOChalky
mol color ColorID 2 ;# grey
mol addrep top

# now the beads will be colored accordingly

set text "name C1A C1B C2A C2B C3A C3B C4A C4B"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 20 ;# green3
mol addrep top

set text "name D2A D3A D2B D3B"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 20 ;# green3
mol addrep top

set text "name NC3"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 16 ;# black
mol addrep top

set text "name PO4"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 19 ;# green2
mol addrep top

set text "name GL0 GL1 GL2"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 15 ;# iceblue
mol addrep top



# PROTEIN

set text "name BB"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 0 ;# blue
mol addrep top

mol selection {name regexp "^(SC).*"}
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 27 ;# magenta
mol addrep top


# CHOLESTEROL

set text "resname CHOL"
mol selection $text
mol representation Licorice 1 12 12
mol material AOChalky
mol color ColorID 3 ;# orange
mol addrep top




display ambientocclusion on
display aodirect 0.400000
color Display Background white
display projection Orthographic


