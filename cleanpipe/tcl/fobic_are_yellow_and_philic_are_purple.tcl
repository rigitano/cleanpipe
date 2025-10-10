#delete all current representations
set num_reps [molinfo top get numreps]
for {set i 0} {$i < $num_reps} {incr i} {
    mol delrep 0 top
}
    
    
    # Create representation for hydrophobic residues (yelow licorice)
    mol addrep top
    mol modselect 0 top {resname ALA VAL LEU ILE MET PHE TRP PRO}
    mol modstyle 0 top Licorice 0.3 12 12
    mol modcolor 0 top ColorID 4
    mol modmaterial 0 top Goodsell

    # Create representation for hydrophilic residues (purple licorice)
    mol addrep top
    mol modselect 1 top {resname ARG ASN ASP GLN GLU HIS LYS SER THR TYR CYS}
    mol modstyle 1 top Licorice 0.3 12 12
    mol modcolor 1 top ColorID 11
    mol modmaterial 1 top Goodsell

    # Add cartoon representation for the backbone
    mol addrep top
    mol modselect 2 top {protein}
    mol modstyle 2 top NewCartoon
    mol modcolor 2 top ColorID 2
    mol modmaterial 2 top Goodsell
