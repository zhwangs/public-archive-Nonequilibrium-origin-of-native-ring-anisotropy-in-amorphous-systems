#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <parent_directory> <start_index> <end_index>"
    exit 1
fi

# always run this for correct output format.

# Assign parent directory and start/end indices from command line arguments
parent_dir="$1"
start=$2
end=$3

# Check if the parent directory exists
# if [ ! -d "$parent_dir" ]; then
#     echo "Error: Directory '$parent_dir' does not exist."
#     exit 1
# fi
if [ ! -d "$parent_dir/out" ]; then
    mkdir "$parent_dir/out"
fi

file_path_out="$parent_dir/out"
python3 pair_distribution.py "$file_path_out"


 
# Loop from start to end index
for ((i=start; i<=end; i+=2))
    
do
    echo "Current '$i' number rings"
    file_path="$parent_dir/n$i"
    file="input_file.dat"

    ct_val="${i}"
    awk -v key="n" -v ct_val="$i" '$1 == key {$2 = ct_val} 1' "$file" > temp && mv temp "$file"

    python3 chain_connect_step1.py "$file_path"
    python3 c_func_step2.py "$file_path"
    python3 cov_closed_loop_step3.py "$file_path"
    python3 convert_full_coord.py "$file_path"

done

python3 ring_distribution_within.py "$parent_dir"
