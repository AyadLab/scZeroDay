#!/bin/bash

source /home/matt/.bashrc

cd /data/mrd/cognitive.seeds_zero/R

if command -v Rscript > /dev/null; then

    echo "Running your R script in the background ..."
    nohup Rscript 03_states.1_2024.05.02.R > 03_script.output_2024.05.02.log 2>&1 &

    # display the process ID of the job
    echo "R script is running with process ID: $!"

else
    echo "Rscript command not found. Please install R or Rscript."
    exit 1
fi

