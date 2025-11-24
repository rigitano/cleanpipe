#!/bin/bash

# Usage example
# ./collect_tprs_within_lambdas.sh nameOfTheFolder

# It will get all the tprs inside the 4_PROD folders within all Lambda folders. 
# That saves all of them in a new folder "collected_tprs" in the current location. 
# Please notice that the prod and bench tprs will be collected.

# Check if the input argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <path_to_search>"
    exit 1
fi

# Directory to search
search_dir=$1

# Destination directory for the .tpr files
dest_dir="collected_tprs"

# Create the destination directory if it doesn't exist
mkdir -p "$dest_dir"

# Temporary file to count the number of collected tprs
temp_count_file=$(mktemp)

# Find and copy the .tpr files
find "$search_dir" -maxdepth 6 -type f -path "*/Lambda_[0-9]*/4_PROD/*.tpr" ! -name '#*' -print0 |
while IFS= read -r -d '' filepath; do

    # Extract the filename without extension
    filename=$(basename "$filepath" .tpr)
    
    # Extract the grandparent directory name and get only the numeric part
    grandparent_dir=$(basename "$(dirname "$(dirname "$filepath")")")
    
   
    new_filename="${filename}.tpr"
  
    # Handle potential filename collisions by appending a counter
    counter=1
    while [ -e "${dest_dir}/${new_filename}" ]; do
        new_filename="${filename}_${counter}.tpr"
        ((counter++))
    done

    # Copy the file to the new directory with the new filename
    echo "collecting ${filepath} to ${dest_dir}/${new_filename}"
    cp "$filepath" "${dest_dir}/${new_filename}"
    echo 1 >> "$temp_count_file"
done

count_tprs=$(wc -l < "$temp_count_file")
rm "$temp_count_file"

echo "${count_tprs} files collected to ${dest_dir}/"

