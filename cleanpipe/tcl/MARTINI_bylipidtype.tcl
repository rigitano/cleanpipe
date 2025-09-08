


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






mol representation Licorice 2.6 50 50
set text "resname CHOL"
mol selection $text
mol material AOChalky
mol color ColorID 11 ;# purple
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname POPC}
mol material AOChalky
mol color ColorID 10 ;# cyan
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DPPC}
mol material AOChalky
mol color ColorID 3 ;# orange
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DSPC}
mol material AOChalky
mol color ColorID 9 ;# pink
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DOPC}
mol material AOChalky
mol color ColorID 7 ;# green
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DOPE}
mol material AOChalky
mol color ColorID 4 ;# yellow
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DOPS}
mol material AOChalky
mol color ColorID 11 ;# purple
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname POPE}
mol material AOChalky
mol color ColorID 15 ;# lightblue
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname POPS}
mol material AOChalky
mol color ColorID 13 ;# lightpink
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname SSM}
mol material AOChalky
mol color ColorID 12 ;# limegreen
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DBPC}
mol material AOChalky
mol color ColorID 14 ;# 
mol addrep top

mol representation Licorice 2.4 50 50
mol selection {resname DLPC}
mol material AOChalky
mol color ColorID 21 ;# 
mol addrep top










# Water – transparent beads
# mol representation VDW 1.0
# mol selection {resname W}
# mol material GlassBubble
# mol color ColorID 0 ;# blue
# mol addrep top

# NA+
# mol representation VDW 1.0
# mol selection {resname NA}
# mol material AOChalky
# mol color ColorID 1 ;# red
# mol addrep top

# CL-
# mol representation VDW 1.0
# mol selection {resname CL}
# mol material AOChalky
# mol color ColorID 7 ;# green
# mol addrep top


puts "Martini lipid visualization applied"