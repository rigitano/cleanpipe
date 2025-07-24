#----------------------------------------------------------------------
# Procedure: get_residue_info
#
# Purpose:
#   For a given molecule (molID), a keyword, and a residue name (rname),
#   this procedure creates an atom selection using the expression:
#       "resname <rname> and <keyword>"
#
#   It then returns a two-element list:
#       {rep_atom_count occurrence_count}
#
#   Where:
#     - rep_atom_count: the number of atoms in one representative occurrence
#       (i.e. for the first distinct residue found).
#     - occurrence_count: the number of distinct residues (using "resid").
#----------------------------------------------------------------------
proc get_residue_info {molID keyword rname} {
    # Create the selection for the given residue name and keyword.
    set sel [atomselect $molID "resname $rname and $keyword"]
    
    # Get unique residue IDs.
    set resids [lsort -unique [$sel get resid]]
    set occurrence_count [llength $resids]
    
    if {$occurrence_count > 0} {
        # Use the first occurrence to get a representative atom count.
        set first_resid [lindex $resids 0]
        set sel_first [atomselect $molID "resname $rname and resid $first_resid and $keyword"]
        set rep_atom_count [$sel_first num]
        $sel_first delete
    } else {
        set rep_atom_count 0
    }
    
    $sel delete
    return [list $rep_atom_count $occurrence_count]
}

#----------------------------------------------------------------------
# Procedure: categorize_residues
#
# Purpose:
#   For each loaded molecule, iterate over a fixed set of keywords.
#   For each keyword:
#     - Compute the overall (molecule+keyword) distinct residue count.
#     - Iterate through all unique residue names in the molecule.
#     - For each residue name that produces a non-empty selection under the
#       keyword, compute the representative atom count (from one occurrence)
#       and the number of occurrences (distinct residue IDs).
#
#   The results are printed to the console.
#----------------------------------------------------------------------
proc categorize_residues {} {
    # Hardcoded keywords (note the last one selects residues that are not
    # any of the other categories).
    set keywords {protein lipid lipids nucleic sugar water waters ion ions "(not protein and not lipid and not lipids and not nucleic and not sugar and not water and not waters and not ion and not ions)"}
    
    # Get the list of loaded molecules.
    set molList [molinfo list]
    if {[llength $molList] == 0} {
        puts "No molecules loaded in VMD."
        return
    }



    
    # Process each molecule.
    foreach molID $molList {
        # Retrieve molecule information.
        set filename   [molinfo $molID get filename]
        set molname    [molinfo $molID get name]
        set numFrames  [molinfo $molID get numframes]
        set numAtoms   [molinfo $molID get numatoms]
        
        puts "\n============================================================"
        puts "Molecule ID:       $molID"
        puts "Molecule Name:     $molname"
        puts "Filename:          $filename"
        puts "Number of frames:  $numFrames"
        puts "Number of atoms:   $numAtoms"
        puts "------------------------------------------------------------"
        
        # Build the complete (unique) list of residue names in this molecule.
        set sel_all [atomselect $molID "all"]
        set resnames [lsort -unique [$sel_all get resname]]
        $sel_all delete
        
        # Iterate over each keyword.
        foreach keyword $keywords {
            puts "Keyword: $keyword"
            
            # For the "double combination" (molecule + keyword):
            # Get all atoms matching the keyword.
            set sel_keyword [atomselect $molID "$keyword"]
            # Get the distinct residues (using "resid").
            set overall_resids [lsort -unique [$sel_keyword get resid]]
            set overall_res_count [llength $overall_resids]
            $sel_keyword delete
            
            # Flag to check if any residue matches under this keyword.
            set any_found 0
            # Iterate over each unique residue name.
            foreach rname $resnames {
                # (Skip empty residue names if any.)
                if {[string length $rname] == 0} { continue }
                
                # Create selection for this triple: molecule, keyword, and resname.
                set sel [atomselect $molID "resname $rname and $keyword"]
                if {[$sel num] > 0} {
                    set any_found 1
                    # Get the representative atom count and occurrence count.
                    set info [get_residue_info $molID $keyword $rname]
                    lassign $info rep_atom_count occurrence_count
                    puts "    resname: $rname (atoms per residue : $rep_atom_count) occurrences: $occurrence_count"
                }
                $sel delete
            }
            if {!$any_found} {
                puts "    none"
            }

            
        }
        puts "\n============================================================\n"
    }
}

#----------------------------------------------------------------------
# Execute the categorization.
#----------------------------------------------------------------------
categorize_residues
