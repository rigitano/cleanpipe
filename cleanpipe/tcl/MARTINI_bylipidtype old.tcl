


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
axes location off

# --- COLOR DEFINITIONS ---
proc vmdrestoremycolors {} {
  set colorcmds {
    {color Display {Background} white}
    {color Display {FPS} white}

    # Martini lipid color scheme (from MartiniGlass / M3 lipidome)
    {color Resname {POPC} cyan}
    {color Resname {DPPC} orange}
    {color Resname {DSPC} pink}
    {color Resname {DOPC} green}
    {color Resname {DOPE} yellow}
    {color Resname {DOPS} purple}
    {color Resname {POPE} lightblue}
    {color Resname {POPS} lightpink}
    {color Resname {CHOL} purple}
    {color Resname {SSM} limegreen}

    # Water beads
    {color Resname {W} blue}
    # Optional: ions
    {color Resname {NA+} red}
    {color Resname {CL-} green}
  }
  foreach colcmd $colorcmds {
    catch {eval $colcmd}
  }
}
vmdrestoremycolors

# --- REPRESENTATIONS ---
# Remove default representation
mol delrep 0 top




# Cholesterol
mol representation Licorice 2.6 50 50
mol selection {resname CHOL}
mol material AOChalky
mol color Resname
mol addrep top




# Lipid headgroups (phosphate, charged beads)
mol representation Licorice 2.6 50 50
mol selection {name NC3 PO4 ROH COO}
mol material AOChalky
mol color Resname
mol addrep top

# Lipid tails (apolar beads, generic rule)
mol representation Licorice 2.4 50 50
mol selection {resname POPC DPPC DSPC DOPC DOPE DOPS POPE POPS and not name NC3 PO4 ROH COO}
mol material AOChalky
mol color Resname
mol addrep top




# Water – transparent beads
mol representation VDW 1.0
mol selection {resname W}
mol material GlassBubble
mol color Resname
mol addrep top

puts "Martini lipid visualization applied"