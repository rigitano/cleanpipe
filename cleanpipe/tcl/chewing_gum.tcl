#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}

# Generate the molecular surface
    mol addrep top
    mol modselect 0 top {protein}
    mol modstyle 0 top Surf
    mol modmaterial 0 top Glossy
    mol modcolor 0 top ColorID 23