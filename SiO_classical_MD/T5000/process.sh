#!/bin/bash

# Define variables
BASE_DIR="$(pwd)" 
NUM_DIRS="$1"
FILE_TO_COPY="md.in"
NUM_ROWS=96
DATA_FILE="sio.data"
" "> "$DATA_FILE"

DUMP_FILE="dump.lammpstrj"


# Check if the file to copy exists
if [ ! -f "$FILE_TO_COPY" ]; then
    echo "Error: File '$FILE_TO_COPY' does not exist in the current directory."
    exit 1
fi

 
# Loop through N directories
for i in $(seq 1 $NUM_DIRS); do
    echo "Current $NUM_DIRS directory."
    # Create directory name
    DIR_NAME="${BASE_DIR}/dir_$i"
    
    
    # Create directory if it doesn't exist
    mkdir -p "$DIR_NAME"
    
    # Copy the file into the directory
    cp "$FILE_TO_COPY" "$DIR_NAME/"
 

    # Change to the directory
    cd "$DIR_NAME" || exit
    # sed -i "/#seed/a insert text here" "md.in"

    cp "/root/Desktop/host/SiO_classical_MD/potential/seed.dat"  "seed.dat"

     awk -v line="$i" 'NR==line' "seed.dat" > "temp_seed.dat"
 
 
 
 
 
    lmp_serial < "$FILE_TO_COPY"
    tail -n "$NUM_ROWS" "$DUMP_FILE" >> "${BASE_DIR}/$DATA_FILE"
 



    rm -rf "$DIR_NAME"

    
    # Change back to the original directory
    cd - > /dev/null

    # Run the program
    #cat  "${BASE_DIR}/$DATA_FILE" 
    
done

echo "Process completed for $NUM_DIRS directories."

