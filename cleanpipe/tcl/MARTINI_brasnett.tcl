

#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}



# Martini bead visualization (tailored to your system)

# --- MATERIALS ---
proc vmdrestoremymaterials {} {
  set mlist {AOChalky GlassBubble AOShiny Glossy}
  set mymlist [material list]
  foreach mat $mlist {
    if { [lsearch $mymlist $mat] == -1 } {
      material add $mat
    }
  }

  # AOChalky – pastel look
  material change ambient   AOChalky 0.0
  material change diffuse   AOChalky 0.85
  material change specular  AOChalky 0.0
  material change shininess AOChalky 0.53
  material change opacity   AOChalky 1.0

  # GlassBubble – transparent water
  material change ambient   GlassBubble 0.25
  material change diffuse   GlassBubble 0.34
  material change specular  GlassBubble 1.0
  material change shininess GlassBubble 1.0
  material change opacity   GlassBubble 0.08
  material change transmode GlassBubble 1
}
vmdrestoremymaterials

# --- DISPLAY SETTINGS ---
display shadows on
display ambientocclusion on
display aoambient 0.9
display aodirect 0.4
display projection orthographic
axes location off



# --- HEADGROUP beads (NC3, PO4) ---
mol representation VDW 1.0
mol selection {name NC3 PO4}
mol color ColorID 21   ;# blue
mol material AOChalky
mol addrep top

# --- BACKBONE beads (GL1, GL2) ---
mol representation VDW 1.0
mol selection {name GL1 GL2}
mol color ColorID 17   ;# yellow
mol material AOChalky
mol addrep top

# --- TAIL beads (C1A...C4B) ---
mol representation VDW 1.0
mol selection {name C1A C1B C2A C2B C3A C3B C4A C4B}
mol color ColorID 7    ;# green
mol material AOChalky
mol addrep top

# --- WATER beads (W) ---
color change rgb 23 0.0 0.5 1.0   ;# define ColorID 23 as blue
mol representation VDW 0.9
mol selection {name W}
mol color ColorID 23
mol material GlassBubble
mol addrep top



# --- CONTECTIONS ---
mol representation Licorice 0.3 12 12
mol selection {not name W}
mol color ColorID 2    ;# grey
mol material AOChalky
mol addrep top




puts "Martini bead visualization applied"
