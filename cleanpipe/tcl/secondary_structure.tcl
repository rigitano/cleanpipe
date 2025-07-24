# Delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}




# Colors for each type of secondary structure


# Beta strands (Extended beta)
mol color ColorID 4 ;# yellow
mol representation NewCartoon
mol selection {backbone and extended_beta}
mol material Opaque
mol addrep top

# Beta strands (Bridge beta)
mol representation NewCartoon
mol color ColorID 17 ;# pale yellow
mol material Opaque
mol selection {backbone and bridge_beta}; # Single hydrogen-bonded beta-strand pairs
mol addrep top

# Alpha helices
mol representation NewCartoon
mol selection {backbone and alpha_helix}
mol color ColorID 1 ;# red
mol material Opaque
mol addrep top

# 3-10 helices
mol representation NewCartoon
mol selection {backbone and helix_3_10}
mol color ColorID 9 ;# pink
mol material Opaque
mol addrep top


# Pi helices
mol representation NewCartoon
mol selection {backbone and pi_helix}
mol color ColorID 13 ;# mauve
mol material Opaque
mol addrep top


# Coils (Unstructured random regions)
mol representation NewCartoon
mol selection {backbone and coil}
mol color ColorID 2 ;# gray
mol material Opaque
mol addrep top

# Turns ( Change direction by a well-defined pattern of H bonds)
mol representation NewCartoon
mol selection {backbone and turn}
mol color ColorID 11 ;# purple
mol material Opaque
mol addrep top










# Now highlight some important amino acids

# Prolines as orange
mol representation Licorice
mol selection {resname PRO}
mol color ColorID 3 ;# orange
mol material Opaque
mol addrep top


# Cysteines (backbone and side chains)
mol representation Licorice
mol selection {resname CYS}
mol color ColorID 2 ;# gray
mol material Opaque
mol addrep top

mol representation VDW
mol color Name
mol selection {element S}
mol material Opaque
mol addrep top




