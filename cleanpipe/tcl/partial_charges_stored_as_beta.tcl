#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}

display projection orthographic


mol selection {all}
mol representation CPK
mol color Beta 
mol material Opaque
mol addrep top