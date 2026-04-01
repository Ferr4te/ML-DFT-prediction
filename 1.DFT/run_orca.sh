#!/bin/bash

ORCA_PATH="/share1/orca/orca"
# For each .inp file
for input in *_gn.inp; do
    # Get the filename without .inp extension
    basename="${input%.inp}"
    
    # Check if output file already exists
    if [ ! -f "${basename}.out" ]; then
        # Run orca
        echo "Starting job for $input"
        $ORCA_PATH "$input" > "${basename}.out" 2>&1
        
        # Wait 1 minute before starting next job
        echo "Waiting 1 minute before next job..."
        sleep 60
    fi
done

echo "All jobs completed!"