


#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}




# --- MATERIALS ---
proc vmdrestoremymaterials {} {
  set mlist {Opaque Transparent Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Goodsell AOShiny AOChalky AOEdgy BlownGlass GlassBubble}
  set mymlist [material list]
  foreach mat $mlist {
    if { [lsearch $mymlist $mat] == -1 } {
      material add $mat
    }
  }

  # AOChalky – pastel/plastic look
  material change ambient   AOChalky 0.0
  material change diffuse   AOChalky 0.85
  material change specular  AOChalky 0.0
  material change shininess AOChalky 0.53
  material change opacity   AOChalky 1.0

  # GlassBubble – for water (transparent beads)
  material change ambient   GlassBubble 0.25
  material change diffuse   GlassBubble 0.34
  material change specular  GlassBubble 1.0
  material change shininess GlassBubble 1.0
  material change opacity   GlassBubble 0.04
  material change transmode GlassBubble 1
}
vmdrestoremymaterials

# --- DISPLAY SETTINGS ---
display shadows on
display ambientocclusion on
display aoambient 0.9
display aodirect 0.4
display projection orthographic
display depthcue off
# axes location off



color Display Background white




# put a grey representation in all, that will be overiden my each specific type
mol representation Licorice 2.3 50 50
set text "not W"
mol selection $text
mol material AOChalky
mol color ColorID 2 ;# grey
mol addrep top


# CHOL
mol representation Licorice 2.6 50 50
set text "resname CHOL"
mol selection $text
mol material AOChalky
mol color ColorID 11 ;# purple
mol addrep top

# PC

mol representation Licorice 2.4 50 50
mol selection {resname POPC}
mol material AOChalky
mol color ColorID 23 ;# blue 2
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DPPC}
mol material AOChalky
mol color ColorID 0 ;# blue
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DSPC}
mol material AOChalky
mol color ColorID 6 ;# silver
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DOPC}
mol material AOChalky
mol color ColorID 15 ;# iceblue
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DBPC}
mol material AOChalky
mol color ColorID 21 ;# cyan 2
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DLPC}
mol material AOChalky
mol color ColorID 10 ;#  cyan
mol addrep top

# SSM

mol representation Licorice 2.4 50 50
mol selection {resname SSM}
mol material AOChalky
mol color ColorID 7 ;# green
mol addrep top


# PE

mol representation Licorice 2.4 50 50
mol selection {resname DOPE}
mol material AOChalky
mol color ColorID 4 ;# yellow
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname POPE}
mol material AOChalky
mol color ColorID 14 ;# ochre
mol addrep top


# PS

mol representation Licorice 2.4 50 50
mol selection {resname DOPS}
mol material AOChalky
mol color ColorID 1 ;# red
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname POPS}
mol material AOChalky
mol color ColorID 9 ;# pink
mol addrep top




# PI

mol representation Licorice 2.4 50 50
mol selection {resname SAPI}
mol material AOChalky
mol color ColorID 3 ;# orange
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname POPI}
mol material AOChalky
mol color ColorID 31 ;# orange 2
mol addrep top





 


#SOLUTION

# Water – transparent beads
# mol representation VDW 1.0
# mol selection {name W}
# mol material GlassBubble
# mol color ColorID 0 ;# blue
# mol addrep top

# NA
# mol representation VDW 1.0
# mol selection {name NA}
# mol material AOChalky
# mol color ColorID 1 ;# red
# mol addrep top

# CL
# mol representation VDW 1.0
# mol selection {name CL}
# mol material AOChalky
# mol color ColorID 7 ;# green
# mol addrep top

# generic name ION
# mol selection "name ION"
# mol representation VDW 1.0
# mol color ColorID 5   ;# tan
# mol material AOChalky
# mol addrep top



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




puts "Martini lipid visualization applied"