# Delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}







#UNSELECTED IS GRAY

# See protein as licorice
mol selection {all}
mol representation Licorice
mol color ColorID 2 ;# grey
mol material Transparent
mol addrep top

# Add a cartoon representation for the backbone
mol selection {protein}
mol representation NewCartoon
mol color ColorID 2 ;# grey
mol material Transparent
mol addrep top



#SELECTED IS RED

# The selection
set text "(resname PRO CYS) or (name CA) or (element S) or (index 1 2 3)"


# represent selection as red Licorice
mol selection $text
mol representation Licorice
mol color ColorID 1 ;# red
mol material Opaque
mol addrep top





# Write index file, to use the selection in gromacs. xxx I dont know if index starts at 0, if they are the same, or if formatting is necessary inside the ndx file

set sel [atomselect top $text]
set ch [open "~/Desktop/custom_selection.ndx" w]
puts $ch [$sel get index]
close $ch



# Clean up
$sel delete
$ch delete
