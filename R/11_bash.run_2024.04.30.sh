#!/bin/bash

source /home/mrd152/.bashrc
cd /data/mrd/spatial/R

# Create output directory if it doesn't exist
output_dir="/data/mrd/spatial/Output"
mkdir -p "$output_dir"

output_dir_log="/data/mrd/spatial/Output/logs/05_deconv_zero"
mkdir -p "$output_dir_log"

fig_dir="/data/mrd/spatial/Output/Figures/05_deconv_zero"
mkdir -p "$fig_dir"

data_dir="/data/mrd/spatial/Output/Rdata/05_deconv_zero"
mkdir -p "$data_dir"

# Loop through each Seurat object
for i in {1..18}; do
  export i="$i" 
  nohup Rscript 05_suter.deconvolution_2024.04.30.R "$i" > "$output_dir_log/output_log_$i.txt" 2>&1 # &
  # wait $! # wait is causing it to not run in bkgd 
done