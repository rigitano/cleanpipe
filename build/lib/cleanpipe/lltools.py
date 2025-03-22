import pandas as pd
import copy
import re

from cleanpipe import bricksStorage


def procv_ll(data_ll, data_key_spec, lookup_ll, lookup_key_spec, lookup_return_cols):
    """
    Performs a "vlookup" join between two tables provided as list-of-lists.
    
    Parameters:
      - data_ll: left table (list of lists)
      - data_key_spec: list of key specifications for the left table.
                       Each element is either an int (column index) or a set of ints.
      - lookup_ll: lookup table (list of lists)
      - lookup_key_spec: key specification for the lookup table.
                         (Same rules as data_key_spec; note that here the numbers refer
                          to the original lookup table column indices.)
      - lookup_return_cols: list of column ids (from the lookup table) that should be returned.
    
    The function has several algorithms:

    #first algorithm:
      1. Converts the inputs into DataFrames.
      2. For the lookup table, renames its columns by prefixing with "L_" so that
         we can later select the desired lookup columns unambiguously.
      3. Constructs a new column "_key" in each DataFrame by applying the key_spec.
         When an element in the spec is a set, a sorted tuple is built.
      4. Drops duplicate keys in the lookup table (keeping only the first match).
      5. Performs an inner merge on the "_key" column.
      6. Returns only the original left table columns plus the selected lookup columns.
         Rows that have no match in the lookup table are omitted.
      7. unmatched lines are stored in unmatched_rows_1 for usage of the next algorithm
    #second algorithm:
      look for cases containing the gromacs anoying X wildcard on the lookupt table






    Examples:
            procv_ll(ll_bonds, [0, 3, 2], ll_bondtypes, [9, 4, 2], [4])
               
               
    
      procv_ll(ll_bonds, [set([3,4]), 2],
               ll_bondtypes, [set([0,1]), 2],
               list(range(3, 10)))
    """


    #translate names to the second algorith to use after the first is done
    ll_boring_table = data_ll
    key_spec_boring = data_key_spec
    ll_lookup_table = lookup_ll
    key_spec_lookup = lookup_key_spec
    l_ids_return    = lookup_return_cols


    ##########check if column ids are out of range#############
    
    def check_key_spec(table, table_name, key_spec):
        """Check that every row in table has enough columns for the given key_spec."""
        for row_idx, row in enumerate(table):
            for spec in key_spec:
                if isinstance(spec, int):
                    if spec < 0 or spec >= len(row):
                        raise IndexError(
                            f"Key spec error in {table_name} at row {row_idx}: index {spec} is out of range "
                            f"(row length: {len(row)}). Row: {row}"
                        )
                elif isinstance(spec, (list, tuple, set, frozenset)):
                    for idx in spec:
                        if not isinstance(idx, int):
                            raise ValueError(
                                f"Key spec error in {table_name} at row {row_idx}: expected integer indices, got {idx} in {spec}"
                            )
                        if idx < 0 or idx >= len(row):
                            raise IndexError(
                                f"Key spec error in {table_name} at row {row_idx}: index {idx} from spec {spec} is out of range "
                                f"(row length: {len(row)}). Row: {row}"
                            )
                else:
                    raise ValueError(f"Unsupported key specification element: {spec}")

    def check_return_indices(table, table_name, return_indices):
        """Check that every row in table has enough columns for the return indices."""
        for row_idx, row in enumerate(table):
            for idx in return_indices:
                if not isinstance(idx, int):
                    raise ValueError(
                        f"Return indices error in {table_name} at row {row_idx}: expected integer index, got {idx}"
                    )
                if idx < 0 or idx >= len(row):
                    return 'problem'
        return 'ok'
                    


    
    # Check indices in both tables before processing.
    check_key_spec(ll_lookup_table, "lookup table", key_spec_lookup)
    check_key_spec(ll_boring_table, "boring table", key_spec_boring)
    result =check_return_indices(ll_lookup_table, "lookup table", l_ids_return)
    if result == 'problem':
        print(f'CLEAN PIPE atention: lookup table begining with {ll_lookup_table[0]} dont have the same number of columns in all lines')




    ######### now the real function ###############

    def construct_key_as_string(row, key_spec):
        """Builds a composite key from a row using the provided key_spec."""
        key_parts = []
        for spec in key_spec:
            if isinstance(spec, int):
                key_parts.append(row[spec])
            elif isinstance(spec, (set, frozenset)):
                key_parts.append(tuple(sorted(row[i] for i in spec)))#sets are converted into sorted tuples. that does the trick
            else:
                raise ValueError(f"Unsupported key specification element: {spec}")
        return str(key_parts)# the key will actually be a string. this is necessary to convert the character X into a regex wilcard



    ######## first algorithm using join ##########
    
    # Convert the left table to a DataFrame; leave column names as their default integers.
    df_data = pd.DataFrame(data_ll)
    
    # Convert the lookup table to a DataFrame and rename its columns (to avoid collisions)
    global df_lookup
    
    df_lookup = pd.DataFrame(lookup_ll)
    n_lookup_cols = df_lookup.shape[1]
    lookup_col_names = ['L_' + str(i) for i in range(n_lookup_cols)]
    df_lookup.columns = lookup_col_names
    
    # Build the key for each row in the left table
    df_data['_key'] = [construct_key_as_string(row, data_key_spec) for _, row in df_data.iterrows()]
    # Build the key for each row in the lookup table using the transformed key spec
    df_lookup['_key'] = [construct_key_as_string(row, lookup_key_spec) for _, row in df_lookup.iterrows()]

    # Drop duplicates in the lookup table (keeping the first match for each key)
    df_lookup_unique = df_lookup.drop_duplicates(subset=['_key'])
   
    # PERFORM THE INNER JOIN (_key)
    merged_df = pd.merge(df_data, df_lookup_unique, on='_key', how='inner')
    
    # Determine the final output:
    #  - All original left table columns (their names are the default integers)
    #  - Plus the lookup columns requested, which we map to the renamed columns.
    data_cols = list(df_data.columns.drop('_key'))
    lookup_return_col_names = ['L_' + str(x) for x in lookup_return_cols]
    
    # Build the final DataFrame and convert it back to a list of lists
    final_df = merged_df[data_cols + lookup_return_col_names]

    # this is the output of the first algothim
    global joined_left_ll
    joined_left_ll = final_df.values.tolist()

    global unmatched_rows_1
    # Compute the rows from the original data table that did not match the lookup (i.e. not present in merged_df)
    unmatched_rows_1 = df_data[~df_data['_key'].isin(merged_df['_key'])].drop(columns=['_key']).values.tolist()




    ######## second algorithm using regex ##########


    
    # Work on a copy
    global output_table
    output_table = copy.deepcopy(unmatched_rows_1)

    # d_lookup = {} #old version using dict
    

    #filter lookup table to keep only rows with the 'X' wildcard
    global ll_lookup_table_only_X
    ll_lookup_table_only_X = [sublist for sublist in ll_lookup_table if 'X' in sublist]

    # loop to add the lookup table to a dataframe. each row of the df will be like '['CA','CB','CG']' | ['0.00','180']. 
    #please notice that the first colum will be a key that is a list and ordered tuples (but converted into a string) and a list with the infos that should be added to the boring table. 
    global df_lookup2
    df_lookup2 = pd.DataFrame(columns=['key','pattern','value'])
    
    for row_idx, row in enumerate(ll_lookup_table_only_X):
        try:
            #this creates the key. something like 'CA', 'CB', ('CD','CE')
            key = construct_key_as_string(row, key_spec_lookup)
        except IndexError as e:
            raise IndexError(f"Error building key for lookup table at row {row_idx}: {e}")
        try:
            return_values = [row[i] for i in l_ids_return]
        except IndexError as e:
            #print(f"less parameters found at row {row_idx}: {e}")
            return_values = [row[i] for i in list(range(l_ids_return[0],len(row)))]

        #now I have to convert the key into a regex pattern. for example '['CA','X']' will become '\[\'CA\',\'[A-Z0-9]+\',\'CG\'\]'
        parts = key.split('X')
        escaped_parts = [re.escape(part) for part in parts] # Escape each literal part so special regex characters are treated literally
        # Join the escaped parts with the wildcard regex pattern
        wildcard = '[A-Z0-9]+'
        pattern = wildcard.join(escaped_parts)

        #this contructs the table
        #d_lookup[key] = return_values #old version with dict
        df_lookup2.loc[len(df_lookup2)] = key, pattern, return_values

    
    # go throught the output table (a copy of the boring table) and for each line, append values found in the lookup table
    for row_idx, row in enumerate(output_table):
        try:
            borring_key = construct_key_as_string(row, key_spec_boring)
        except IndexError as e:
            raise IndexError(f"Error building key for boring table at row {row_idx}: {e}")
            
        #this line is the essence of all. this is the lookup. using wildcards

        # make a reverse regex lookup of sorts. It will see if any patterm in the column pattern in the lookup dataframe matches the current key
        result = df_lookup2.loc[df_lookup2['pattern'].apply(lambda pat: re.fullmatch(pat, borring_key) is not None), 'value']
     
        # extract a value from the 'value' column for the first matching row:
        if not result.empty:
            found_value = result.iloc[0]
        else:
            found_value = ["Unknown"]

        #put the result in the output list of lists
        row.extend(found_value)



    
    
    return joined_left_ll + output_table #concatenate results from the first algorithm to the results of the second



def clean_comments_out(list_of_lists):
    """
    Given a list of lists (rows), this function removes any columns
    after a cell containing a semicolon (";").
    
    If a cell (which is a string) contains a semicolon mixed with text,
    the cell is truncated to only the text before the semicolon (with
    surrounding whitespace removed). All subsequent cells in that row are dropped.
    Cells that are not strings (like numbers) are kept as-is.
    
    Parameters:
        list_of_lists (list of list): The parsed data.
        
    Returns:
        list of list: The cleaned data.


    example
    cl.clean_comments_out(ll_bonds)
    """
    cleaned_rows = []
    for row in list_of_lists:
        new_row = []
        for cell in row:
            if isinstance(cell, str) and ';' in cell:
                # Truncate cell at the first semicolon and remove extra spaces.
                new_cell = cell.split(';', 1)[0].strip()
                if new_cell not in [";",""]:
                    new_row.append(new_cell)
                # Stop processing further columns in this row.
                break
            else:
                new_row.append(cell)
        cleaned_rows.append(new_row)
    return cleaned_rows


def split_ll_diherals_into_proper_and_improper(ll_dihedrals):
    """
    give a list of lists representing dihedrals, this function will split them into two list of lists (proper and improper)
    
    example usage
    ll_dihedrals_proper , ll_dihedrals_improper           = scl.plit_ll_diherals_into_proper_and_improper(ll_dihedrals)
    """


    # First list: sublists where the 4th element equals '2'
    ll_improper = [sublist for sublist in ll_dihedrals if len(sublist) > 3 and sublist[4] == '2']
    # Second list: sublists where the 4th element is not '2'
    ll_proper = [sublist for sublist in ll_dihedrals if len(sublist) > 3 and sublist[4] != '2']
    return ll_proper, ll_improper




def filter_dihedrals_to_keep_only_phi_and_psi(ll_dihedrals,ll_atoms):
    """
    filter to keep only backbone dihedrals
    
    given ll_dihedrals, will find out the ones that the ones that define phi and psi
    to find that, ll_atoms is necessary, because the only way to find out is to check if the connection between N-CA (using C as 0) ad CA-C (using N as 0)

    example usage:

    ll_atoms     = dll_parsed_molecules['Protein_chain_A']['[ atoms ]']
    ll_dihedrals = dll_parsed_molecules['Protein_chain_A']['[ dihedrals ]']
    
    ll_backbone_dihedrals = cl.filter_dihedrals_to_keep_only_phi_and_psi(ll_dihedrals,ll_atoms)
    
    """


    df_hihedrals = bricksStorage.ll2df(ll_dihedrals)



    df_atoms = bricksStorage.ll2df(ll_atoms)
    df_atoms = df_atoms.iloc[:, :5] # keep only 'id','ff name','resid','resname','structural name'
    df_atoms.columns=['id','ff name','resid','resname','structural name']
    
    
    filtered_rows = []
    for index, row in df_hihedrals.iterrows():
        # Check if the atoms involved in the dihedarl are in phi and psi aroung the CA. this are the possibilities: 
        # ["C N CA C","C N CA C","N CA C N","N C CA N"]
        #if they are, we add the line to the filtered table. if they are not, the loop goes on without adding the line
        #all the nested ifs are to avoid unecessary lookups
        name_i = df_atoms.loc[df_atoms['id'] == row[0], 'structural name'].values[0]
        if name_i == "C":
            name_l = df_atoms.loc[df_atoms['id'] == row[3], 'structural name'].values[0]
            if name_l == "C":
                name_j = df_atoms.loc[df_atoms['id'] == row[1], 'structural name'].values[0]
                name_k = df_atoms.loc[df_atoms['id'] == row[2], 'structural name'].values[0]
                if (name_j == "N" and name_k == "CA") or (name_j == "CA" and name_k == "N"):
                    #add
                    filtered_rows.append(row)
                    continue
        else:
            if name_i == "N":
                name_l = df_atoms.loc[df_atoms['id'] == row[3], 'structural name'].values[0]
                if name_l == "N":
                    name_j = df_atoms.loc[df_atoms['id'] == row[1], 'structural name'].values[0]
                    name_k = df_atoms.loc[df_atoms['id'] == row[2], 'structural name'].values[0]
                    if (name_j == "C" and name_k == "CA") or (name_j == "CA" and name_k == "C"):
                        #add
                        filtered_rows.append(row)
                        continue
    
    
    df_filtered = pd.DataFrame(filtered_rows)
    
    
    return bricksStorage.df2ll(df_filtered)